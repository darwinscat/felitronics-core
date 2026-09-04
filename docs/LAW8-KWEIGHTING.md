<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# Law 8 — silence was costing 54× more than music

[`DSP-ARCHITECTURE.md`](DSP-ARCHITECTURE.md) law 8 says every feedback kernel flushes denormal state in
**software**, because hardware FTZ is a desktop optimisation and the `wasm-audio` tier has no FP control
register at all — a browser cannot set it. Two kernels were not obeying it. This is the measurement, the
mechanism, and why the fix is placed where it is.

---

## The defect

`analysis::KWeightingFilter` — the BS.1770 pre-filter, two TDF-II biquads in `double` — had no flush. With
zero input its state decays into the subnormal range and then **sticks**: for the 48 kHz RLB stage
(`a1 = −1.99004745483397971`, `a2 = 0.99007225036620994`) the pair `(z₁, z₂) = (−251u, +249u)`, with
`u = 2⁻¹⁰⁷⁴`, maps to itself **exactly** —

| check | value | rounds to | needed |
|---|---|---|---|
| `a₂·M` | 248.508135 | 249 | 249 ✓ |
| `a₁·M + N` | −250.501911 | −251 | −251 ✓ |

— so the state never reaches zero on its own. Measured out to 600 s of silence, frozen at `−1.240105e-321`,
bit-identical on arm64 and x86.

**The cost, measured twice on different benchmarks.** x86-64, i9-13900H, gcc 14.2 `-O2 -ffp-contract=off`,
2 s of signal, best of 5:

| | music (440 Hz) | silence (stuck state) |
|---|---|---|
| no flush | 0.299 ms | **16.077 ms** — 54× |
| no flush, hardware FTZ/DAZ forced on | 0.304 ms | 0.304 ms — no penalty at all |
| per-block flush | 0.299 ms | 0.384 ms |

The middle row is the point: forcing FTZ on the *same binary* removes the gap entirely, so this is purely
the denormal path. On Apple Silicon, which runs subnormals at full speed, both columns read the same and
nothing is visible — which is why it survived.

The suite's own gate now measures it as a ratio: on x86 the fix takes silence/music from **74.85× to
1.00×**. On arm64 it reads ~0.9× either way and proves nothing there, which is stated in the test.

---

## Why the flush is per **sub-hop**, and not where it first went

The first implementation flushed at the end of `LoudnessMeter::process()`. That is wrong twice over, and
both reasons were found by the review consilium rather than by writing it.

**1. It would break chunk invariance.** The end of `process()` is wherever the *caller* chose to cut the
stream. Putting a numerical event there makes the answer depend on the host's block size — and this repo
explicitly claims and tests the opposite, bit-exact across call sizes 1 … 100 003
(`tools/tests/ProbeTests.cpp`). The 10 ms sub-hop is a deterministic amount of *audio*, so the cadence is
identical however the caller chunks.

**2. The binding constant is 90 ms, not 2.85 s — and it is the shelf, not the RLB.** From the 1e-15 flush
threshold, the two stages reach the subnormal floor at very different speeds:

| rate | RLB (pole 0.995) | **shelf (pole 0.856)** | 10 ms sub-hop |
|---|---|---|---|
| 44.1 kHz | 125 650 samples (2.85 s) | **3 983 samples (90 ms)** | 441 samples |
| 48 kHz | 136 778 (2.85 s) | **4 324 (90 ms)** | 480 samples |
| 192 kHz | 548 142 (2.85 s) | **17 342 (90 ms)** | 1 920 samples |

I had computed only the RLB and concluded there was a 15× margin against `core::kMaxBlockSize = 8192`.
That was the wrong stage. 8192 samples is 170 ms at 48 kHz — already **twice** the shelf's 90 ms — so a
per-host-block flush would leave the shelf sitting subnormal through half of every silent block. Against
the sub-hop the margin is 9×, at every rate, by construction.

`LoudnessMeter::process()` also never clamps `n`, and `fcore::Probe` documents rates down to 1 kHz, where a
single 8192-sample call spans 8.2 seconds. A cadence tied to the caller cannot be made safe.

---

## What it changes, stated honestly

**Not bit-exact, and the claim was withdrawn rather than defended.** On the CI parity fixture the *reported*
LUFS and dBTP are bit-identical before and after (`-0x1.4275e021a3d7cp-2`, `0x1.71011e4ec516p+1`), and 5 of
97 pre-gate block energies change — all in the decaying silence tail, between **−458 dB and −1290 dB**,
against an absolute gate at −70 LUFS. Four of the five become exactly zero.

But one fixture is not a theorem. A short gap followed by programme is a constructed counter-example: the
flush deletes a residual δ ≈ 1e-15 that would otherwise still be ringing when the music returns, and the
cross-term survives into block energies **above** the gate at a relative size around 1e-15 (≈ 1e-14 dB).
So the defensible statement is:

- reported loudness moves by less than 1e-13 dB;
- pre-gate block energies agree to ~1e-12 relative;
- native↔wasm parity stays **bit-exact**, because both tiers run the same flush at the same deterministic
  boundaries — which the CI NULL test re-confirms.

---

## The second kernel

`multiband::MultibandProcessor` never called `splitter_.flushDenormals()`. The method has existed on
`eq::MultibandSplitter` since it was written; nothing in `modules/multiband/` invoked it, so the crossover
SVFs decayed into the subnormals on silence and stayed. The per-band `Compressor` flushes itself, which is
precisely why this went unnoticed. One line, after the split loop, plus a test whose negative control fails
without it.

---

## Still open, deliberately

- **A single non-finite sample still freezes the meter forever.** K-weighting is an IIR, so its state stays
  NaN; `NaN > absT` is false, the absolute gate silently drops every later block, and the meter keeps
  reporting a healthy number computed over the pre-NaN part while `droppedBlocks()` says nothing. Using
  `flushPoison` here would be free and would end the freeze — but a meter should *report* damage, not
  quietly repair it, so the real fix needs an observability counter (`nonFiniteSubHops()`) and a rewrite of
  the tests that currently pin the freeze. That is its own change, not a rider on a performance fix.
- **`fcore::Probe` accepts sample rates the loudness path cannot serve.** Its floor is 1 kHz, but the shelf
  is designed at 1682 Hz, so below **3364 Hz** the bilinear `tan(π f₀/fs)` goes negative and the shelf pole
  leaves the unit circle — |pole| = 1.31 at 3 kHz, **2.15 at 1 kHz**. The filter is unstable there and the
  loudness answer is meaningless, while `prepare()` happily accepts it.
- **Other unflushed recursive state**, found in the same audit and left alone on purpose:
  `saturation::Saturator`'s DC blocker flushes non-finite only, *by design*, with a written rationale that
  deserves its owner's decision rather than a drive-by change; `core::Smoother` can stick when its target is
  exactly zero (a mute), though it has no in-repo caller of `next()`; `poweramp`'s coefficient smoothers and
  `rigplayer`'s gain state are the same shape. The claim that K-weighting was "the only feedback kernel
  without a flush" was true only of `core` + `analysis` + `oversampling`, and is not repeated here.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# Analyzer performance — where a multi-resolution tick goes, and `MultiResSpectrumPaneFast`

**Status:** ✅ landed · **Module:** `felitronics::analysis`
**Sibling docs:** [`ANALYZER-MULTIRES.md`](ANALYZER-MULTIRES.md) is the physics. This is the cost.

---

## TL;DR

A `MultiResSpectrumPane` tick was **44 % transcendental math**, and most of it was work nobody asked
for: a `log10` for every bin whose only consumer was an `exp` that undid it, `exp2` on a loop
invariant, the display tilt computed twice per column, and the entire geometry of a column derived
once for the fill trace and again for the peak.

`MultiResSpectrumPaneFast` removes it. **Same numbers, no approximation** — the fill is bit-identical
and NULL-tested as such. It is a **sibling**: `MultiResSpectrumPane` keeps its behaviour and its
tests, and the two are compared against each other in
`modules/analysis/tests/MultiResFastPaneTests.cpp`.

| | Apple M5 Pro | i9-13900H |
|---|---|---|
| current pane, scalar FFT | 402.0 µs | 740.1 µs |
| **fast pane, scalar FFT** | **335.6 µs** (1.20×) | **553.8 µs** (1.34×) |
| current pane, pffft | 149.9 µs | 324.3 µs |
| **fast pane, pffft** | **82.2 µs** (1.82×) | **139.3 µs** (2.33×) |

The saving is the **same number of microseconds on both backends** — 66 µs on the M5, 185 µs on the
i9 — because what it removes does not depend on the transform. That is the whole thesis, and it is
why the ratio looks unimpressive on the scalar FFT (which is dominated by the transform) and large on
a SIMD one (which is not). It is also why the saving is ~2.8× bigger on x86: glibc's `log10`, `exp`
and `pow` cost far more than Apple's.

---

## The profile that found it

`perf record` on the i9, a driver doing nothing but ticks of the pane on pffft, `ingest` and
`buildColumns` in separate `noinline` phases:

```
54.79 %  the pane itself (incl. 6.08 % of harness noise generation)
43.98 %  libm.so.6            <-- forty-four percent of runtime is transcendental math
 1.08 %  libc.so.6
```

Attributed by phase and symbol:

| phase | share | of which libm |
|---|---|---|
| `ingest` | 62.5 % | `log10` 20.6 %, `exp` 4.9 % |
| `buildColumns` | 26.3 % | `log10` 2.1 %, `pow` ≈2.3 %, `exp2` 2.8 % |

The band integration everyone assumed was the cost is **7.6 %** (`readPower` + `bandPower`). The FFT
itself is **4.6 %** (`radf4_ps`). Neither is the problem.

A crew seat instrumented the real header and counted the calls per tick at the default ladder:

| phase | `exp2` | `log2` | `pow` | `log10` | `exp` |
|---|---|---|---|---|---|
| `ingest` | 0 | 0 | 0 | **10 755** | **10 755** |
| `buildColumns` | **6 348** | 2 986 | 2 703 | 3 564 | 0 |
| `starve` | 0 | 0 | 0 | 0 | 10 755 |

**37 111 calls per tick, of which 1 762 survive the changes below** — and those 1 762 are the one
honest `log10` per emitted value.

## What changed, and why each is exact

1. **The peak trace lives in power.** The pane stored the hold in dB, which required `10·log10(pw)`
   for every bin every tick to produce the value the peak law compares against, and then
   `exp(pk·ln10/10)` for every bin to convert the hold back to power for its prefix sum — a transform
   and its own inverse, 10 755 of each. In power, "fall `peakFallDb` per tick" is a multiply by the
   constant `10^(−peakFallDb/10)`. `max()` commutes with a monotone map, so this is the same trace.
   The per-bin dB array goes with them; `tierBinDb` computes its one value on demand.
2. **One column plan, shared by both traces.** `readDb` and `readPeakDb` differ by a single pointer,
   yet each re-derived the whole geometry of the column: the owning tier, the seam, the blend width
   and weight, the band's fractional bin edges, the display tilt. None of it depends on the frame. It
   is computed once per column into a fixed-size plan and applied to both prefix arrays; a band
   becomes four prefix loads and a multiply-add. The plan is rebuilt only when its inputs move — a
   resize, a new sample rate, a tuning edit, a tier becoming valid — never per tick. The bin-limited
   branch is evaluated in the sibling's own operation order (`p / span * cell`, not `p * (cell/span)`)
   so the fill stays bit-identical.
3. **DC and Nyquist are peeled** out of the magnitude loop, leaving a contiguous interior. Bin for
   bin the arithmetic is unchanged.

## What was considered and rejected

A lossy variant was designed in parallel and **not built**. Three independent reviews reached the same
conclusion: after the exact work above, the classic lossy move — a bit-trick `log2` for the remaining
per-column `toDb` — is worth about 1.5 % and would become the largest term in the NULL tolerance. The
measured rejections, kept here so they are not re-proposed:

- **float32 prefix sums.** Not a saving; a dynamic-range loss. A band is a *difference* of two partial
  sums, so beside a full-scale neighbour the float32 cancellation floor sits near −69 dBFS: a −80 dB
  band reads the floor. Worse, a noise-only NULL passes, because noise has no loud partial sum to
  cancel against. The doubles are load-bearing.
- **Column decimation / a coarser log-f grid.** Measured 4–40 dB errors on bin-centred tones, and it
  saves nothing in the case that matters (eight panes at 360 columns, where a 1/24-octave band already
  spans 1.4 columns).
- **Whole-bin band edges.** A 3 dB sawtooth at exactly the seams, which the design note already
  predicts.
- **Deriving the peak from the fill.** `Σ max` is not `max Σ`; measured 1.8–3.4 dB apart on a ramp.
  A different law, not an approximation of this one.
- **Multirate lower tiers** (decimate before the FFT — the 16384 tier only serves up to ~1.02 kHz, so
  ÷16 gives the same 341 ms aperture and the same 2.93 Hz bin from a 1024-point transform) is the one
  lossy idea with a large win left. It needs ~−107 dB stopband rejection to keep a −90 dB probe inside
  budget beside a 0 dBFS blocker. Not attempted here.

---

## The tables

One UI tick = `ingest` + the columns. 48 kHz, frame 16384, hop 1600, default ladder 16384/4096/1024.
Every variant is measured in **one process, round-robin** — each round times each row once and the
figure is the median over nine rounds — so the boost/thermal drift of a laptop is charged to every row
equally. **Absolute microseconds depend on the machine's thermal state; the ratios between rows do
not**, and the ratios are what this table is for. Reproduce with
`felitronics_spectrum_pane_perf_tests` (scalar) and `felitronics_pffft_tests` (both backends), which
print the same pair and fail if the fast pane ever stops being cheaper.

Machines: **Apple M5 Pro** (macOS, AppleClang, `-O2`, NEON) and **i9-13900H** (Debian 13, gcc 14.2,
`-O2`, SSE, `performance` governor, pinned to a P-core with `taskset`).

### 900 columns — one large pane

| pane | M5 Pro µs | i9 µs | M5 ×8 %core | i9 ×8 %core |
|---|---|---|---|---|
| classic 1024, scalar | 20.4 | 47.2 | 0.49 | 1.13 |
| classic 2048, scalar | 34.8 | 72.9 | 0.83 | 1.75 |
| classic 4096, scalar | 65.6 | 123.8 | 1.57 | 2.97 |
| classic 8192, scalar | 129.0 | 240.3 | 3.10 | 5.77 |
| classic 16384, scalar | 262.2 | 489.5 | 6.29 | 11.75 |
| multi-res current, scalar | 402.0 | 740.1 | 9.65 | 17.76 |
| **multi-res FAST, scalar** | **335.6** | **553.8** | **8.06** | **13.29** |
| classic 1024, pffft | 11.5 | 34.4 | 0.28 | 0.83 |
| classic 2048, pffft | 15.1 | 43.4 | 0.36 | 1.04 |
| classic 4096, pffft | 23.4 | 60.1 | 0.56 | 1.44 |
| classic 8192, pffft | 39.5 | 100.9 | 0.95 | 2.42 |
| classic 16384, pffft | 74.4 | 180.1 | 1.79 | 4.32 |
| multi-res current, pffft | 149.9 | 324.3 | 3.60 | 7.78 |
| **multi-res FAST, pffft** | **82.2** | **139.3** | **1.97** | **3.34** |

On x86 the fast multi-resolution pane (139.3 µs) is now **cheaper than a single classic 16384 pane**
(180.1 µs) — three tiers of constant-Q for less than one flat-resolution transform. The classic pane
has not had this treatment and still pays the per-bin `log10` and the per-column `pow`.

### 360 columns — the eight-pane case (an amp view)

`%core` is `8 × tick / 33333 µs`: eight panes all ticking in the same 30 fps frame on one thread.

| pane | M5 %core | i9 %core |
|---|---|---|
| classic 2048, pffft | 0.25 | 0.65 |
| classic 16384, pffft | 1.65 | 3.80 |
| multi-res current, pffft | 2.91 | 6.07 |
| **multi-res FAST, pffft** | **1.87** | **2.91** |
| multi-res current, scalar | 8.94 | 16.14 |
| **multi-res FAST, scalar** | **8.03** | **12.95** |

Without pffft a full-ladder multi-resolution analyzer is still not affordable for eight panes on x86.
Shipping the SIMD backend is worth more than any of this, and a shorter tier ladder
(`setTiers`) is worth more than both.

---

## The denormal bug this work uncovered (fixed separately, v0.23.x)

On an all-zero frame the fill smoother is `pw ← (1−c)·pw`, and in float32 that descent does not end at
zero: at the smallest subnormal `c·pw` rounds away and the state **sticks**, leaving every bin doing
subnormal arithmetic on the message thread. Nothing sets FTZ/DAZ there, so on x86 the pane got
**slower the longer the transport stayed stopped** — 505.4 µs on settled silence against 167.6 µs with
FTZ forced, a **3.0× penalty**, while on ARM it was free and therefore invisible. Fixed by flushing the
smoothed power to a true zero below `1e-30` — deliberately **not** `core::flushDenormal`, whose `1e-15`
is an amplitude threshold and would erase the −150…−200 dB bins this pane deliberately sums.
The classic `SpectrumPane` never had it: it smooths dB, which is already floored.

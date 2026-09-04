<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# One bad sample used to make the loudness meter lie quietly

`analysis::LoudnessMeter` survived a non-finite input sample by reporting a healthy, plausible number for a
programme it had stopped measuring. This is the fix, and the two things found while making it.

---

## The defect

K-weighting is an IIR. One NaN sample made its state NaN **forever**, so every later 400 ms block energy was
NaN — and `NaN > absT` is **false**, so the absolute gate silently dropped all of them. The meter then
averaged only the part of the programme that predated the NaN and answered with it. `droppedBlocks()` said
nothing: it counts capacity overflow.

Measured on a 5 s 1 kHz tone with one NaN at 1 s:

| | before | after |
|---|---|---|
| non-finite block energies | **40 of 47** | 0 |
| `integratedLufs()` | −12.034538 — plausible, and wrong | −12.034502 |
| `momentaryLufs()` | **−120 forever** (reads as "bypassed") | recovers within 400 ms |
| the event is reported | nowhere | `nonFiniteSubHops() == 1` |

The error before was not merely larger — it was **unbounded**, because it depended on where the NaN fell.
Early in the programme the meter reported −120 LUFS, silence, for loud material. After the fix the damage is
one 10 ms sub-hop recorded as silence: 1/40 of a gating block, `10·log10(39/40) = −0.11 dB` on the blocks
that overlap it, and a few thousandths of a dB on the integrated answer. Bounded and derivable.

## The fix, in three parts

1. **Heal.** `KWeightingFilter::flushDenormals()` now uses `core::flushPoison` — the strict superset that
   also zeroes non-finite state — at the same deterministic 10 ms boundary law 8 already flushes on.
2. **Catch the energy, per channel.** The poisoned channel's sub-hop is recorded as silence. Per *channel*
   and not on the weighted sum, because BS.1770 gives LFE `w = 0` and `0 · NaN` is NaN: one excluded channel
   would otherwise poison a sub-hop whose audible channels were fine.
3. **Report.** `nonFiniteSubHops()`, sticky until `reset()`, beside `droppedBlocks()`. Non-zero means the
   number above is best effort, not a measurement.

Also: `setChannelWeight()` now refuses a non-finite weight. Healing the filter cannot heal a poisoned
*configuration* — `w[c] · energy` would make every sub-hop non-finite from a caller's mistake.

**Why the sub-hop is the unit.** It is the first energy boundary that does not depend on how the caller
chunks its calls, and it is the quantum momentary, short-term and LRA already share. A 400 ms block counter
was the other candidate and is worse twice over: one bad sub-hop lands in one to four overlapping blocks
depending on phase, so a single event reports as 1, 2, 3 or 4; and a programme shorter than 400 ms never
forms a block at all, so a 300 ms clip containing a NaN would report zero.

## Two things found on the way

**A `+inf` sample reached `(int)` in the LRA histogram — undefined behaviour, on `main`.** `NaN > absT` is
false, but `+inf > absT` is **true**. An infinite short-term energy passed the absolute gate, made the
relative threshold infinite, passed `inf >= inf`, and computed `(int) ((lufsOf(inf) + 70) · 10)`. Reproduced
with UBSan: *"inf is outside the range of representable values of type 'int'"*, from one `+inf` sample
aligned so its sub-hop was the last in a 3 s window. The clamp two lines later is too late — the conversion
itself is the UB. Every gate now tests `isfinite` first, and the case is pinned in the suite that the
sanitizer job runs.

Healing made this *more* reachable, not less: with the state repaired the neighbouring energies are finite,
so a single infinity can sit alone in a window instead of being swamped by NaNs. Hardening the gates was
therefore not optional.

**The energies are recorded as `0.0`, not kept raw — and that is the opposite of what it looks like.**
Keeping the poison reads as the honest choice, since `gatingBlockEnergies()` is the forensic surface. But it
is *also* the `wasm-audio` tier's cross-toolchain bit-comparison surface, and a **computed** NaN is not
portable. `inf − inf` — which is exactly what a `+inf` input produces inside the biquad — is

| | arm64 (Apple clang) | x86-64 (gcc 14.2) |
|---|---|---|
| `inf − inf` | `0x7ff8000000000000` | `0xfff8000000000000` |
| a NaN *literal* | `0x7ff8000000000000` | `0x7ff8000000000000` |

measured on both machines. On `origin/main`, **37 of 41 block-energy lines differ between the two
architectures** for one `+inf` input; with this change they are identical. A NaN that *arrived* as input
carries its payload identically on both — which is precisely why the old pinned test never showed it.

## The test that pinned the defect could not see it

The group in `tools/tests/ProbeTests.cpp` existed to make this behaviour impossible to change unnoticed. It
did not notice: its 8192-sample lead-in is 171 ms, shorter than one 400 ms gating block, so `lufsBefore` was
the −120 no-block sentinel and the "frozen" assertion compared a sentinel with itself. It passed both before
and after the behaviour changed. It is rewritten around a fixture long enough to have a reading, and the
`+inf` group's "reads as SILENCE" assertion is gone — that signal now reads −29.29 LUFS.

## Verified

macOS/Apple clang 76/76 · gcc 14.2 on Debian 13 76/76 · `wasm-audio` 75/75 in node · ASan+UBSan clean ·
native↔wasm parity bit-identical · **healthy audio bit-identical to `origin/main`** on the parity fixture,
so nothing about clean streams changed.

## Still open

Layer 1 of the finding: rejecting non-finite input at the **C-ABI boundary**, which belongs to the facade
(P9). It does not replace this: a desktop plugin never crosses that boundary — a NaN arrives from another
plugin upstream — and an RT plugin cannot refuse to emit a buffer.

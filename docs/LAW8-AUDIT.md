<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# Law 8 — the whole audit, closed

[`DSP-ARCHITECTURE.md`](DSP-ARCHITECTURE.md) law 8: *every feedback kernel flushes denormal state in
SOFTWARE; hardware FTZ is a desktop optimisation, never a correctness crutch.* Three passes have now
touched this. The first fixed `analysis::KWeightingFilter` and `multiband::MultibandProcessor`
([`LAW8-KWEIGHTING.md`](LAW8-KWEIGHTING.md)) and left a short "still open, deliberately" list. This
document is the pass that finishes it: **every recursive kernel in the repo has a verdict here, including
the ones that turned out clean and why.** Nobody should have to derive any of it a fourth time.

The short version: **six kernels were stalling, four are clean, and one of the four "leftovers" the first
pass recorded as harmless was in fact the most expensive of the lot.**

---

## The mechanism, stated generally (it is not really about denormals)

A one-pole approaching a target never *arrives*: `x ← t + r·(x−t)`. Under round-to-nearest the residual
has a **band of fixed points** — write the residual as `k` ulps and it maps to itself whenever

```
k·(1 − r) ≤ ½          i.e.   k ≤ ½ / (1 − r)
```

so the decay stops dead as soon as it enters that band. Two consequences, and the audit needed both:

- **Target 0.** The band is in the *subnormal* range, where the ulp is absolute (2⁻¹⁰⁷⁴ double, 2⁻¹⁴⁹
  float). The state parks there **permanently** and every later operation touching it pays the subnormal
  penalty — 10–100× on x86 without FTZ, nothing at all on Apple Silicon, which is why this class of bug
  survives on the dev machine. This is law 8 proper.
- **Any other target.** The band is at `k·ulp(target)`, which is a *normal* number, so there is no CPU
  penalty — but the kernel still never reports itself settled, and anything gated on exact settling
  (`core::Smoother`'s cost-zero early-out) never fires. Same defect, invisible cost.

A flush written only against the first case fixes half the problem. Both fixed kernels that glide toward a
*parameter* — `core::Smoother`, `poweramp`'s coefficient smoothers, `rigplayer`'s gain ramps — therefore
carry **two** conditions: a magnitude threshold, and "the update did not move the value", which is exact,
target-independent, and cannot fire early (in the normal range the mantissa makes `k(1−r) ≫ ½`).

---

## Verdicts

| kernel | recursive state decaying to 0? | sticks? | cost of not fixing | verdict |
|---|---|---|---|---|
| `core::Smoother` | yes, asymptotic to target | **yes** — 720 ulp (3.56e-321) after 22 s of `next()`; 5 ulp after 28 s of `advance(128)` | one subnormal multiply **per sample**; plus a `std::pow` every block, forever, on a settled parameter | **FIXED** |
| `core::LinearSmoother` | no — fixed-increment ramp | — | — | clean by construction |
| `measurement::LevelProbe` | no persistent state | — | — | clean |
| `measurement::Sweep` | no recursion at all | — | — | clean |
| `dither::Dither` | loop exists, state does not decay | no | — | clean |
| `saturation::Saturator` DC blocker | yes | **yes** — 2.14e-42 from 1.5 s, permanent | **23.9× measured**, and it leaks out of the stage | **FIXED** |
| `rigplayer` blend / trim ramps | yes | **yes** — 2 ulp, stall from ~1 s | ~3 assisting ops **per sample per channel**, forever | **FIXED** |
| `poweramp` coefficient smoothers | yes | **yes** — ~5 ulp | `topoCur`/`leakCur` reach the **per-oversampled-sample** kernel | **FIXED** |
| `analysis::CorrelationMeter` | yes | **yes** — 7200 ulp after ~220 s | three subnormal double FMAs **per sample** | **FIXED** |
| `dynamics::AutoLeveler` | yes | **yes** — 28 ulp after ~106 s | two FMAs **per block**, behind a silence gate | **FIXED** (see "not independently tested") |

### The six that were stalling

**`core::Smoother`** — the one the previous pass named. `current = target + coeff·(current−target)`; with
`target == 0` (a mute, a band parked at 0 dB) it is a pure geometric decay. At 30 ms / 48 kHz,
`coeff = 0.9993058`, so every `k ≤ 720` is a fixed point: measured stuck at **720 ulp after 22.17 s**,
first subnormal at 21.25 s. `advance(n)` sticks whenever `k(1−coeffⁿ) ≤ ½` — 5 ulp at n = 128, and it does
reach zero at n = 1024 only because `coeff¹⁰²⁴ = 0.49 < ½`. In-repo caller: `eq::EqBand`, 15 instances,
`advance()` only.

The second half of this one is not a denormal at all. The header claimed "an idle parameter costs one
compare per block, not a transcendental (the eq lanes' cost-zero contract rides on this)", guarded by
`fabs(current − target) > 0.0`. A residual that merely decays never becomes `0.0`, so the guard never
fired and the `pow` ran forever. A 1e-15 absolute snap alone does **not** fix that for a real target:
measured with threshold-only, target 1000 at n = 64 still made **18750 `pow` calls over 25 s** with a
residual of 1.25e-12. The fixed-point condition is what makes the claim true, for every target.

**`saturation::Saturator`'s DC blocker — the one the first pass got wrong, including this author.** It
flushed non-finite state only, with a written rationale: *"no denormal threshold here: zapping a tiny
finite tail would break the bit-identical NULL for quiet legitimate signals."* Three things turned out to
be true at once:

- The stall is **permanent**, not transient. It was assumed the state would traverse the subnormals in
  ~0.25 s and underflow to zero. It does not: `R = 0.99967283` (10 Hz at 4× 48 kHz) gives a fixed-point
  band up to `k ≤ 1528`, and `y1` freezes at **2.14e-42 after 1.49 s** and stays there. On silence the
  Asym curve emits exactly zero (`tanh(k·bias) − biasTanh_` is the same float expression), so nothing
  re-excites it.
- **It escapes the stage.** The frozen value feeds the 128-tap downsampling FIR with every operand
  subnormal, and the Saturator's own output becomes **−3.40249e-41 on every sample** — bit-identical on
  arm64 and x86 — i.e. contagion into whatever comes next.
- **Measured cost**, x86-64 (i9-13900H, gcc 14.2 `-O2`, 2 s mono, best of 5): music **18.97 ms**, stalled
  silence **452.72 ms — 23.9×**, 22.6 %RT per channel. Forcing FTZ/DAZ on the *same binary*: 18.33 ms, no
  penalty at all, which is what proves it is purely the denormal path. With the flush: 18.25 ms.

The rationale is retracted rather than edited, because the NULL it appealed to never tested the claim
either way: scenario 4's 1e-20 block is absorbed to exact zero by the Asym bias (`−0.3f + 1e-20f == −0.3f`;
`ulp(0.3) = 3e-8`), so the fixture never puts a tiny *finite* value into that state. Verified rather
than assumed: three threshold variants were built against it — per-sample 1e-30f, per-sample 1e-15f,
per-call 1e-15f — and all **197** checks of the three saturation suites pass in every one, so the
fixture does not discriminate the threshold at all. (194 before this document's own law-8 test added
three; if you reproduce this, count what your tree has rather than trusting the figure.) The threshold
chosen is **1e-30f**,
on the margin argument alone: 8 decades above the float subnormal floor, ~600 dB below anything audible.
The flush is applied **per sample, not per call** — a threshold at the end of a call would put a numerical
event wherever the *caller* cut the stream, which is exactly the chunk-invariance argument from
[`LAW8-KWEIGHTING.md`](LAW8-KWEIGHTING.md), and this module asserts bit-identical chunking explicitly.

**`rigplayer`'s blend and trim ramps.** `end = want + (current−want)·decay`, per block. An exact-zero
target needs no exotic pack: `BlendKnob::linOf` returns `0.0` for any level at or below −120 dB, which is
how a pack spells "this path is off" — the repo's own test rig ships it (dry end `wetDb −120`, wet end
`dryDb −120`). Move the dial there after an audible position and the gain decays to a subnormal fixed
point (`decay = 0.766` at block 128 → `k ≤ 2.1`, so 2 ulp; stall from ~1 s) and stays, while the mix loop
— gated on `dryActive_`, a dry IR being **loaded**, not on either gain — keeps running
`gd += stepDry; a[c][i]·gw + d[c][i]·gd` over it, ~3 assisting operations per sample per channel, for the
life of the rig. What does *not* stall is the initialised or reset `0.0f`: `0 → 0` stays exactly 0.

**`poweramp`'s 13 block-rate coefficient smoothers.** Most are consumed behind a per-block gate
(`presOn`, `depthOn`, `loadOn`, `biasOn`, `ironOn`, all `> 1e-4f`), which reads a stuck subnormal as
"off" — the right answer, reached by accident. **`topoCur` is not gated**: `TubeStage` blends
`(1−topo)·pp + topo·se` on every *oversampled* sample, so after one SE→PP toggle a stuck topo costs a
subnormal multiply and add per sample per channel at 4× rate, forever; `leakCur` reaches the per-sample
curves the same way. So this was never the 13-FMAs-per-block housekeeping it looked like. Snapping also
restores a real property: `topo` lands on exact 0, so the "all-off ⇒ bare push-pull path, byte-identical"
contract holds after a toggle and not only from a cold start. The per-*sample* states in this file
(`dcx1/dcy1`, `otLp/otHf`) were already flushed at 1e-30f and are clean.

**`analysis::CorrelationMeter`** — found in this pass, and it is the first pass's own lesson repeating:
`flushDenormals()` existed and **nothing in the repo called it**. It is a leaf primitive with no in-repo
owner (only a comment in `stereo::StereoWidth` mentions the class), so the adapter driving it was the only
thing that could have, and never did. Three double one-poles decay on silence and stick (`k ≤ 7200` at the
default 300 ms window). The stall is slow — ~220 s of silence to reach the subnormal band — but that is a
long tail on a mix bus, not an exemption. `process()` now flushes itself, so there is no call to forget;
the cost is three `fabs`+compares against a three-FMA dependency chain the loop is latency-bound on anyway.

**`dynamics::AutoLeveler`** — two per-block mean-square followers, stick at `k ≤ 28` after ~106 s. Two
FMAs per block behind a silence gate, so the cost is genuinely immeasurable; fixed because the remedy is
two compares and law 8 is not a cost-benefit rule.

### The four that are clean, and why

- **`core::LinearSmoother`** — a fixed-increment ramp *arrives*. `countdown` hits 0 and the final step
  assigns `currentValue = target` exactly. No asymptote, no residual, nothing to flush. This is the
  structural difference from the exponential `Smoother` beside it.
- **`measurement::LevelProbe`** — the one-pole is a *function-local* `lp`, re-zeroed on every `fill()`
  call and driven by full-scale LCG noise the whole time; there is no persistent state and no silence.
  Off-thread/offline besides.
- **`measurement::Sweep`** — clean twice over: closed-form `sin(φ)·env` with no recursion of any kind, and
  explicitly offline / message-thread (double, allocates).
- **`dither::Dither`** — the error-feedback loop is real, but the state cannot decay. `e[]` is in LSB
  units and is re-excited every sample by a ±1 LSB TPDF draw: with `in == 0`,
  `e = code + fb` where `code` is a small integer and `|fb| ≤ Σ|h_k|·8` (16 for the Weighted curve, 171
  for the 9-tap), so `e` is O(1). Two independent reasons beyond that: for the Weighted curve the
  coefficients 1.5 / −0.5 are **dyadic**, so every `e` is a dyadic rational with a bounded denominator and
  *cannot* be subnormal; and auto-blanking clears the shaper to exact zero after 4096 samples of digital
  black. The `autoBlank == false` path was checked separately and is clean for the same reason.

### Swept and clean, no action (recorded so the sweep is not repeated)

`TruePeakLimiter` (its `grDb` is flushed) · `dynamics::RelativeLevel` (dB domain, fixed point is nonzero)
· `NoiseGate` (per-sample flush plus an env flush) · `eq::EqBand` (per-block `flushState()`, and a reset
when no lane runs) · `eq::Svf` / `eq::MatchedBiquad` (1e-15f) · `rigplayer::runBands` (flushes each
`MatchedBiquad` per block) · `nam::BlendLaw` (linear step + clamp lands on exact 0/1) · the convolver IR
crossfade (a countdown, not a decay) · `stereo::MonoBass` (flushes its crossover, resets on bypass) ·
`stereo::StereoWidth` (`LinearSmoother` only) · `MultiResSpectrumPane` (has its own power-aware guard) ·
`SagEnvelope` (1e-30f) · `dynamics::EnvelopeFollower` in **Rms** mode (1e-30f — the state there is a
POWER, so the house 1e-15 amplitude threshold would have zapped a level of 3.2e-8, i.e. −150 dBFS. That
is not merely audible-adjacent, it made the OUTPUT depend on the caller's block partition, because the
flush fires once per `process()` call: measured on a −160 dBFS input, −7.5 dB of gain reduction in one
10000-sample call against 0.00 dB in 10000 one-sample calls, the same stream. Same reasoning, and the
same constant, as `SagEnvelope` and the power-aware guard in `MultiResSpectrumPane`) · `TubeStage` (memoryless) · `poweramp`'s `gApplied`/`postApplied` (linear ramps,
land exactly) · every `EnvelopeFollower` owner (Compressor, NoiseGate, DeEsser, DynamicEqBand,
LaneDynamics, TransientShaper) · `DelayLine` / `DryAligner` / `StreamResampler` (no recursion) ·
`Fft` / `Pffft*` / `MatrixConvolver*` / `IrResampler` / `PolyphaseOversampler` / `BlendKernels` (FIR and
transforms — no decaying recursion) · `TruePeakMeter` and `SpectrumPane` (settled in the F1 pass) ·
`measurement` / `blend` / `lineareq` / `io` (offline or FIR).

**Not verified:** NAM/neural inference state (external code; an LSTM's zero-input fixed point is nonzero,
so no subnormal stall is expected, but it has not been measured).

### ⚠ THE SWEEP ABOVE ASSUMED A SMALL BLOCK, and one call shape falls outside it

Every "clean, per-block flush" verdict here rests on `core/FlushToZero.h`'s own reasoning: *"a block is
too short to re-traverse the gap, so state never reaches subnormal across blocks."* **That is false for a
big block**, and a big block is not exotic — `Compressor`, `EqEngine` and `Dither` all accept any length,
and `TruePeakLimiter`'s header explicitly tells an offline caller that sizing `maxBlock` to a whole file
is a normal thing to do. A flush that fires once per CALL cannot fire inside one.

Measured (P6, `eq::EqEngine`, one +6 dB bell, 1000 samples of tone then 39000 of digital silence, tail
`[2000, 40000)`):

| call shape | exact zeros | SUBNORMAL | last non-zero |
|---|---|---|---|
| whole file, one call | 1 | **37678** | 7.006e-45 at i = 39999 |
| block 4096 | 35905 | 1774 | 7.006e-45 at i = 4095 |
| block 64 | 38000 | 0 | — |
| block 1 | 38000 | 0 | — |

So the entire tail sits in the subnormal range — the 10–100× stall this law exists to prevent, and the
one already measured at 23.9× on `Saturator`'s DC blocker. `eq::EqBand` and `stereo::MonoBass` are listed
as clean above; that verdict holds for a host-sized block and NOT for a whole-file call. `Saturator` is
genuinely clean either way, because its flush is per SAMPLE.

**What to do about it is a choice, and P6 took the third option:** (a) move the flushes into the sample
loop — correct, but a versioned behaviour change in a module TabbyEQ ships; (b) cap the block a stage
will accept — reintroduces the silent-truncation defect P2 F2 closed; (c) **stop the caller's block from
reaching the stage at all**, which is what `felitronics::mastering` does with a fixed internal quantum,
and why it has one. A future re-audit should measure the whole-file shape for every row above rather than
assume the sweep covered it.

---

## How these are tested

State tests, in the style of `modules/dynamics/tests/DynamicsTests.cpp:61` — after a long silence the
state is **exactly** zero. Not timing: PR #107 removed CPU-grading from the suite and none of it comes
back here. The 23.9× above is a one-off measurement for this write-up, run by hand on x86, not a ctest.

Every one of these was verified by **reverting its own fix and watching the assertion fail** — the flush
is never called from the test, so a test that passes with the flush removed is not testing it.

| kernel | observable | negative control (fix reverted) |
|---|---|---|
| `core::Smoother` | `value()` after 5 s aimed at 0, via `next()` and `advance(128)`; plus a **non-zero** target reached exactly and `settled(0)` firing | 4.15e-73 — a normal double, so FTZ cannot mask it |
| `analysis::CorrelationMeter` | after 5 s of silence, 1000 samples on L only ⇒ `sRR == 0` ⇒ `d == 0` ⇒ `correlation()` is exactly 1.0 | residual 1.9e-19 gives `d = 1.6e-10`, **161× the meter's own 1e-12 gate**, and `correlation()` reads ~1e-9 |
| `saturation::Saturator` | the **output** — a state that reached zero emits exact zero | −3.4e-41 on every sample |
| `rigplayer` | `liveWet()`, the gain the audio thread actually applied | 1.8e-35 at 0.8 s |

Two of them assert **twice**, at an early instant and at 5 s, and the reason is worth keeping: at 5 s the
un-flushed value is *subnormal*, so a machine running hardware FTZ would read it as zero and a regression
could hide. The early assertion is placed where the un-flushed value is still a **normal** float
(Saturator: the flush fires at 1.10 s, the value goes subnormal at 1.39 s, so 1.25 s discriminates on any
machine; rigplayer: flush at 0.69 s, subnormal at 0.87 s, assert at 0.8 s).

**Not independently tested, deliberately:** `dynamics::AutoLeveler` and `poweramp`'s *gated* coefficient
smoothers. Both are invisible through their public surface — the AutoLeveler's silence gate returns before
the state can influence anything, and every poweramp gate reads a stuck subnormal as "off", which is the
same answer as a snapped zero. There is no assertion that would fail without the fix, so writing one would
be theatre. `poweramp`'s `topoCur` is the exception in principle (it is consumed per sample), but at
`u = 0` both `pp` and `se` are exactly zero and with signal the subnormal is absorbed, so it is not
output-observable either.

---

## Checking a new kernel

1. Does it hold state updated as `x ← f(x)` where `f` contracts toward a fixed point? If not (FIR,
   transform, closed-form generator, countdown, linear ramp), stop — it is clean, and say so in the header
   so the question is not re-asked.
2. Can the target be **exactly zero** on a real path? Trace it: a default, a reset, a clamp floor, a dB
   conversion that underflows (`≤ −120 dB` in this repo returns literal `0.0`).
3. Compute the fixed-point band, `k ≤ ½/(1−r)`. If `k ≥ 1` the decay stops; it does not "eventually
   underflow".
4. Ask **where the stuck value is read**, not just where it is written. The three expensive findings here
   were all cheap at the write site and costly at the read site — a per-sample loop, an oversampled
   kernel, a 128-tap FIR.
5. Flush inside the kernel, not in a method the owner must remember to call. That mistake has now been
   made twice in this repo (`MultibandSplitter`, `CorrelationMeter`) and caught twice by an audit rather
   than by a test.
6. Where the state glides to a *parameter* rather than to zero, use both conditions (threshold **and**
   "the update did not move the value"). A threshold alone leaves every non-zero target parked a few ulp
   short forever.

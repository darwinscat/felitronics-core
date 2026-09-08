# What `core::StreamResampler` costs — measured

`felitronics::core::StreamResampler` rate-matches a live stream in and out of a rate-locked neural
model (`nam::NamStage` holds **two** instances per channel — `down` host→model before the network and
`up` model→host after it), and it engages only when `|hostSR − modelRunSR| > 0.5` (`NamStage.cpp`).
**At a 48 kHz host with a 48 kHz capture there is no resampler in the path at all**, so everything
below belongs to 44.1 kHz and other off-native rates.

> ## 🔴 THE KERNEL CHANGED IN P34 — READ THIS FIRST
>
> Up to **v0.26.0** this class was a **Catmull-Rom cubic with no anti-aliasing filter of any kind**.
> P32 measured what that cost; P34 replaced it with a **64-tap polyphase Kaiser-windowed sinc
> (β = 8.6, cutoff 0.99 of the lower Nyquist, 512 phase rows with linear interpolation)**.
>
> **Everything in sections 1–4 describes the CUBIC.** Those measurements are not wrong and they are
> not obsolete — they are the reason the kernel changed, and they are the "before" column of the
> acceptance. They are simply no longer a description of the shipped code. Section 5 is the decision;
> **section 6 is what shipped and what it measures.** Section 7 is the latency consequence, rewritten.
>
> | | cubic (≤ v0.26.0) | 64-tap sinc (P34) |
> |---|---|---|
> | carrier at 17.64 kHz, one round trip | **−4.17 dB** | **+0.0002 dB** |
> | worst phase at 17.64 kHz | **−9.27 dB** | **+0.0000 dB** |
> | modulation depth, 20 kHz | 14.0 dB | 0.0004 dB |
> | decimation stopband, 22.1 → 23.9 kHz | **−3.0 dB flat, 0.0 dB sample peak** | −9.1 … −88.8 dB |
> | a 20 kHz tone's demodulated 100 Hz line, high-gain capture | **−17.65 dBFS (+14.5 dB over its own carrier)** | **−96.54 dBFS (−60.7 dB under it)** |
> | round-trip delay at 44.1 kHz | 3.84 samples (0.087 ms) | **61.40 samples (1.392 ms)** |
> | %RT, one mono round trip, arm64 | 0.0129 | 0.1447 |
>
> **This is an audible change on every model at 44.1 kHz**, in the direction the capture was made in.
> See `CHANGELOG.md` under Unreleased, marked BREAKING (behaviour).

The header used to justify the old kernel with one sentence — *"the driven nonlinear stage masks the
interpolation images"* — and no number. This document is the number, and the number retired the
sentence.

## 1. The mechanism: one thing, not two — *(applies to ANY phase-dependent kernel, then and now)*

A phase-dependent interpolation kernel is a **linear periodically time-varying** filter. For a
complex exponential the output is the ideal output times a per-sample complex gain that depends only
on the interpolation phase `t`:

```
y[k] = e^{jΩ·i_k} · Σ_j w_j(t_k) e^{jΩ·o_j}          o_j ∈ {−1, 0, 1, 2}
     = ideal · M(t_k),   M(t) = Σ_j w_j(t) e^{jΩ(o_j − t)}
```

`M` is 1-periodic in `t`, and its Fourier coefficients are `c_k = H(Ω + 2πk)` — the kernel's own
continuous transform sampled at the image frequencies. So the **amplitude modulation and the
interpolation images are the same phenomenon in two coordinates**; there is no separate "modulation"
defect to argue about, and no separate "images" defence.

Split `M(t) = A₀ + M̃(t)`. `A₀` is the **coherent carrier** — a plain frequency-response error on the
wanted signal. `M̃` is the **added** part: the sidebands, i.e. the images. Masking is a statement
about added components only; it cannot apply to `A₀`.

## 2. THE CUBIC — both axes, round trip 44.1 ↔ 48 kHz *(history: this is what was replaced)*

Coherent carrier / worst phase, dB. Three independent oracles agree to 0.01 dB: an exact
rational-phase enumeration of the closed form, a brute-force simulation, and the shipped class driven
through `NamStage`'s own plumbing (cos and sin runs combined into a per-sample complex gain — the
class is linear, so no bucketing is needed and none is done).

| f, Hz | coherent | worst phase | best phase | **×2: coherent** | **×2: worst** | **×2: best** |
|---|---|---|---|---|---|---|
| 1 000 | −0.00 | −0.00 | −0.00 | −0.00 | −0.00 | −0.00 |
| 5 000 | −0.05 | −0.08 | −0.01 | −0.09 | −0.13 | −0.06 |
| 8 000 | −0.28 | −0.50 | −0.05 | −0.56 | −0.79 | −0.34 |
| 10 000 | −0.64 | −1.16 | −0.11 | −1.30 | −1.82 | −0.77 |
| 12 000 | −1.22 | −2.28 | −0.21 | −2.53 | −3.59 | −1.48 |
| 15 000 | −2.59 | −5.14 | −0.41 | −5.48 | −8.01 | −3.16 |
| **17 640** | **−4.17** | **−9.27** | −0.61 | **−9.03** | **−13.16** | **−5.20** |
| 19 000 | −4.98 | −12.20 | −0.69 | −10.91 | −15.45 | −6.30 |
| 20 000 | −5.48 | −14.79 | −0.74 | −12.09 | −17.85 | −7.02 |

**The ×2 columns are the shipped OrbitCab chain, and they are not the first three doubled.** OrbitCab
runs TWO `NamStage`s in series at the host rate — preamp (`src/core/CabEngine.cpp`) → EQ → poweramp
router (same file) — and `src/PluginProcessor.cpp` says so in its own comment: *"The two NAM stages
each rate-match independently, so their latencies SUM (possible future optimisation: … a single
round-trip instead of two)."* So the shipped path at 44.1 kHz is **four** Catmull-Rom stages, not two.

Doubling the decibels would give −8.35 / −18.53 at 17.64 kHz. The cascade measures **−9.03 / −13.16**:
the carrier is *worse* than doubling and the worst phase *better*, for the same reason the inherited
one-round-trip table was wrong — each stage converts part of the previous one's sidebands back onto the
carrier. What actually changes character is the **best** phase: it falls from −0.61 dB to −5.20, so the
top octave is attenuated at *every* phase instead of only at some, and the modulation depth narrows
from 8.7 dB to 8.0 while the whole band sinks. Two oracles agree on every cell (the shipped class in
`NamStage`'s plumbing, and an independent brute-force double simulation); the composite period is still
147 output samples.

**The set is complete, and it belongs to the shipped priming.** At 44100/48000 = 147/160 the composite
gain is periodic with **exactly 147 output samples**: stage 2's position advances by 160/147 per output,
so after 147 outputs it advances by exactly 160 — a whole number of stage-1 phase periods, since stage 1
advances by 147/160 and repeats every 160. (Checked numerically: deviation from period-147 repetition is
2.6e-13 over five periods, no proper divisor of 147 is a period, and the statistics over 23 520 samples
are identical.)

What those 147 pairs are **not** is the full 160 × 147 product. The stage-1 phase at the point stage 2
reads is a deterministic function of the stage-2 phase, so the cascade traces a **line** through that
torus, and which line is set by how the two stages are primed — here identically, `pos = 1` and `len = 3`,
both reset together in `NamStage::configureRates`. Swept over all 160 integer alignments:

| f | shipped coherent / worst | coherent over all alignments | worst over all alignments |
|---|---|---|---|
| 15 000 | −2.59 / −5.14 | −2.36 … −3.43 | −3.53 … −5.33 |
| 17 640 | −4.17 / −9.27 | −3.59 … −6.83 | −6.96 … −9.29 |
| 20 000 | −5.48 / −14.79 | −4.47 … −13.14 | −13.17 … −15.70 |

So the shipped priming happens to land **at or near the alignment-wide worst** — 0.02 dB off at
17.64 kHz, 0.91 dB off at 20 kHz — which is an observation about this priming, not a bound: the worst
phase itself ranges over 2.3 dB (17.64 kHz) and 2.5 dB (20 kHz) across alignments. The **coherent
carrier is an interference term** — stage 2 converts part of
stage 1's sidebands back onto the carrier, and how much depends on the alignment. The mean over alignments
is exactly the product of the two single-stage carriers (−5.06 dB at 17.64 kHz, −7.77 at 20 kHz), which is
what a "two independent stages" estimate would give and why that estimate is not the cascade.

## 3. THE CUBIC — the decimating direction had no stopband *(history)*

Going 48 → 44.1 the kernel must remove everything above 22.05 kHz, or it folds. It removes
essentially nothing — at phase `t = 0` the Catmull-Rom weights are `(0, 1, 0, 0)`, a bare sample
pick, which attenuates nothing at any frequency:

| tone at 48 kHz | shipped | 16-tap sinc | 32-tap sinc | 64-tap sinc |
|---|---|---|---|---|
| 22 100 Hz | **−3.0 / −0.0** | −10.0 / −9.4 | −15.7 / −15.7 | −34.2 / −34.2 |
| 22 500 Hz | **−3.0 / −0.0** | −11.4 / −10.3 | −19.9 / −19.9 | −53.7 / −53.7 |
| 23 000 Hz | **−3.1 / −0.0** | −13.2 / −11.4 | −26.4 / −26.4 | −88.6 / −85.4 |
| 23 500 Hz | **−3.1 / −0.0** | −14.8 / −12.1 | −34.7 / −34.4 | −93.5 / −91.7 |
| 23 900 Hz | **−3.1 / −0.0** | −15.4 / −12.4 | −42.1 / −39.5 | −92.0 / −90.6 |

dB rms / sample peak, unit input. The two columns measure different things and both are real. No
single spectral LINE exceeds −4.67 dB; the sample peak reaches 0 dB for a time-domain reason — at
`t = 0` the weights are `(0,1,0,0)`, so that output simply IS an input sample, and the phase returns
to within ~1e-13 of 0 every 147 outputs, so one output in every 147 reproduces an input sample to
float rounding (measured output peak 0.999657 against an input peak of exactly 1.0). What survives does not land in one place: a tone at *g* ∈ (22.05, 24) kHz
comes back as **two** strong components — `44100 − g` at about −5 dB and `g − 3900` at about −7 dB —
plus weaker terms near −45 dB. The two STRONG components land in **18.15–22.05 kHz**, not just the top
slice; the weak ones go lower still (14.3 kHz at −42.6 dB for a 22.1 kHz input), so that band is where
the damage is, not where it ends. Measured
on the shipped kernel, 2-second coherent window:

| tone in | components out |
|---|---|
| 22 100 Hz | 22 000 Hz −4.67 · 18 200 Hz −7.92 · 18 100 Hz −47.8 · 14 300 Hz −42.6 |
| 23 000 Hz | 21 100 Hz −5.33 · 19 100 Hz −7.04 · 17 200 Hz −45.9 · 15 200 Hz −43.2 |
| 23 900 Hz | 20 200 Hz −6.06 · 20 000 Hz −6.23 · 16 300 Hz −44.5 · 16 100 Hz −44.2 |

## 4. THE CUBIC — what the driven nonlinear stage actually did to it *(history, and the reason for P34)*

Measured on the reference NAM DI through real captures, against an exact band-limited reference
conversion (`shipped − ideal`, both at 44.1 kHz):

* the **output** leg — the half with no possible causal masker, since nothing after
  `up.produceExact()` in `NamStage` or `RigPlayer` is nonlinear — is the **larger** half, by 5…10 dB
  in every run;
* against the **model's own aliasing floor** (what a WaveNet at 48 kHz folds down by itself, measured
  on an arm with no resampler in it), the **output leg alone** sits below it on a high-gain capture
  (3–24 dB) but above it on a clean one at every level from 17.5 kHz up — **up to +9.8 dB, which is at
  19 kHz and −30 dBFS**. (At 15.25 kHz it is above the floor at two of the three levels, by up to
  +3.8 dB.) Counting the whole rate-match,
  including the model reacting to a droop-ed input, it is above the floor in **22 of 30** tone × level
  cells, by up to **+12.3 dB** — below only on the high-gain capture at 12.3–15.25 kHz;
* driving harder does not help, and in the audible midrange it actively hurts. Per band, error over
  ideal-output across a 42 dB input sweep (high-gain capture, real DI):

  | band | −24 dB drive | 0 dB | +18 dB | change |
  |---|---|---|---|---|
  | 19–22 kHz | −5.3 | −6.1 | −5.8 | flat |
  | 16–19 kHz | −11.0 | −10.9 | −10.5 | flat |
  | 12–16 kHz | −16.8 | −16.0 | −15.3 | +1.5 |
  | 4–8 kHz | −36.5 | −33.5 | −31.3 | +5.2 |
  | 1–4 kHz | −51.7 | −45.2 | −41.2 | **+10.5** |
  | 0–1 kHz | −63.4 | −53.2 | −49.1 | **+14.3** |

  Where the resampler's own artifacts live (16–22 kHz) the ratio is FLAT — masking would require it to
  improve. Where it grows is the bass and low mid, and nothing the resampler does lands there: the
  folded-alias term alone stays below −59 dB in every band under 16 kHz at every drive. What grows is
  the nonlinearity **demodulating** the input leg's images down into the audible range.
* and that demodulation is the loudest single thing measured in this whole task. A 20 kHz tone at
  −18 dBFS into a high-gain capture, shipped path: **a 100 Hz line at −17.7 dBFS — 14.5 dB louder than
  the 20 kHz carrier that produced it** (−32.2). The same bin through the output leg alone is −81.8,
  and through an ideal round trip −174.6. All the IMD products of the image lattice `20000 ± 3900m`
  fall on a 100 Hz grid, which is where they land. The clean capture does it too, 21.4 dB quieter
  (−39.1 dBFS at 100 Hz against the high-gain −17.7). A steady 20 kHz tone is not guitar, so read this as the mechanism rather than
  as a level; the level on real DI is the table above.

So the old one-line defence fails three ways: it never covered the carrier droop at all, it is
conditional where it does apply (3–24 dB below the model's own floor on a driven capture, up to +9.8 dB above
it on a clean one, output leg alone), and the driven stage does not mask the input leg's artifacts —
it moves them into the part of the spectrum where nothing masks anything.

## 5. The decision: which kernel, and why the cutoff and not the tap count decided it

A polyphase windowed-sinc (Kaiser β = 8.6, `L` taps, `P` phases, linear interpolation between
phases) fits the same seam and the same RT contract (table built in `reset()`). Round trip, coherent
/ worst phase. At 32 and 64 taps the two columns **coincide** — the modulation is gone and only a
designed band edge remains; at 16 taps they still differ by 0.04–0.28 dB, i.e. it shortens the
modulation without removing it:

| kernel | 10 k | 15 k | 17.64 k | 19 k | 20 k | %RT (mono, both stages) | round-trip delay |
|---|---|---|---|---|---|---|---|
| shipped Catmull-Rom | −0.64/−1.16 | −2.59/−5.14 | −4.17/−9.27 | −4.98/−12.20 | −5.48/−14.79 | 0.0126 | **3.84 samples** |
| L=16, cutoff 0.94 | −0.00/−0.00 | −0.20/−0.20 | −2.30/−2.34 | −5.30/−5.44 | −8.86/−9.14 | 0.0372 | 15.36 |
| L=32, cutoff 0.94 | −0.00/−0.00 | 0.00/0.00 | −0.11/−0.11 | −1.77/−1.77 | −6.10/−6.10 | 0.0790 | 30.71 |
| **L=64, cutoff 0.99** | **−0.00/−0.00** | **0.00/−0.00** | **0.00/0.00** | **0.00/0.00** | **−0.01/−0.01** | 0.1549 | 61.41 |
| L=64, cutoff 0.99, **×2** | −0.00/−0.00 | 0.00/−0.00 | 0.00/0.00 | 0.00/0.00 | −0.03/−0.03 | 0.3098 | 122.8 |
| L=96, cutoff 0.99 | −0.00/−0.00 | 0.00/−0.00 | −0.00/−0.00 | 0.00/0.00 | 0.00/0.00 | 0.2450 | 92.06 |

UP-leg stopband (rms/peak, dB), since a wider passband buys its transparency with near-Nyquist
rejection: L=32/0.94 gives −15.7 / −19.9 / −26.4 / −34.7 / −42.1 at 22.1 / 22.5 / 23 / 23.5 / 23.9 kHz;
L=64/0.99 gives −9.1 / −15.4 / −27.5 / −47.6 / −88.8; L=96/0.99 gives −11.0 / −22.9 / −52.6 / −88.1 /
−89.9. Against the shipped −3 dB flat, any of them is a different order of thing.

### 🔴 The two tone axes do NOT decide this — program material does

The obvious pick from the table is the cheapest kernel that flattens both axes, i.e. L=32. On real DI
through a driven capture it is **not a clean win**. Rendered through the same model against the same
ideal reference, with the LTI part (the kernel's own band edge — a deliberate design choice) separated
from the residual, error over ideal output per band:

| kernel | drive | 0–1 k | 1–4 k | 4–8 k | 8–12 k | 12–16 k | 16–19 k | 19–22 k |
|---|---|---|---|---|---|---|---|---|
| shipped Catmull | 0 dB | −53.21 | −45.38 | −33.86 | −24.86 | −18.26 | −14.65 | −9.13 |
| L=32, cutoff 0.94 | 0 dB | **−48.47** | **−42.41** | −33.98 | −27.61 | −20.93 | −14.11 | −12.42 |
| **L=64, cutoff 0.99** | 0 dB | **−53.97** | **−48.01** | **−39.63** | **−33.18** | **−26.50** | **−19.34** | **−13.70** |

L=32 at a 0.94 cutoff is **5 dB WORSE in the bass** than the kernel it would replace: its band edge sits
at 20.7 kHz and removes content the reference keeps, and the driven model converts that difference into
low frequencies by exactly the mechanism of §4 — the same demodulation from a different cause. Widen the
cutoff to 0.99 and lengthen the kernel to keep a stopband, and it is better in **every band at every
drive** (checked at −24, 0 and +18 dB).

**So the cutoff is the deciding variable here, not the tap count** — the same axis P31 has open for the
oversampler. A candidate judged on the two tone axes alone picks the wrong one.

**And on the shipped OrbitCab chain the case is stronger than one round trip suggests, in both
directions.** Through two round trips the shipped kernel reaches −9.03 dB coherent and −13.16 worst at
17.64 kHz (−12.09/−17.85 at 20 kHz), while the L=64 / 0.99 candidate stays at 0.00 / −0.03 — it does
not accumulate at all, which is what a flat passband means. But its latency accumulates exactly:
**2 × 61.4 = 122.8 host samples, 2.78 ms** at 44.1 kHz, and its CPU doubles to 0.3098 %RT (≈5 % of what
two `NamStage`s spend on their models). The single-round-trip figures elsewhere in this document
describe `rigplayer`, which runs its two `NamStage`s in PARALLEL, not in series.

The %RT column above is one machine (Apple Silicon, Apple clang); absolute figures do not travel, so
what matters is the fraction. Measured on two, with the same standard WaveNet capture at 44.1 kHz:

| machine | whole `NamStage` (model included) | shipped resampler | 32-tap sinc | what that swap costs |
|---|---|---|---|---|
| arm64, Apple clang | 5.3–6.1 %RT mono | 0.0126 | 0.0790 | **+1.2 %** of the stage |
| x86-64 i9, gcc 14.2 | 11.0–11.6 %RT mono | 0.0515 | 0.1952 | **+1.3 %** of the stage |

The two-axis, stopband and delay numbers are identical on both toolchains to the digits printed here.
The recommended L=64 / 0.99 costs **0.1549 %RT** on the same Mac, i.e. **+2.5 %** of what the stage
already spends on the model — and **61.4 host samples, 1.39 ms**, against today's 3.84 samples
(0.087 ms). That latency is the whole price, and it is the part only Oleh can weigh.

`rigplayer` consequences, checked in the code rather than assumed: slot alignment runs on
`AlignmentTable::delayOf()` → `blendDelay()` / `lagTail_`, and **none of them reads
`latencySamples()`** — `RigPlayer::latencySamples()` only republishes the max of the two slots
outward. So a kernel change moves the absolute PDC and nothing else inside the player, as long as
both slots take the same delay. `AlignmentTable::measureAlignment()` takes the rate as an argument
and defaults to 48 kHz, where the resampler is bypassed entirely; called at another rate it still
measures the same processed samples.

## 6. WHAT SHIPPED — the 64-tap polyphase sinc, measured (P34)

The kernel is `kTaps = 64` taps of a Kaiser-windowed sinc (β = 8.6) sampled into `kPhases = 512` phase
rows built once in `reset()`, with linear interpolation between the two rows bracketing the wanted
phase. The window is a fixed 64 **input** samples wide in both directions; what adapts to the ratio is
the **cutoff**, `fc = 0.99·min(1, outRate/inRate)` in units of the input Nyquist — so both legs of the
44.1 ↔ 48 kHz round trip band-limit to the same **21 829.5 Hz**.

Why that construction and not a time-stretched one: it is the one whose numbers the §5 candidate table
actually describes. Rebuilt from scratch and driven through the same plumbing, it reproduces the §5
up-leg stopband row to better than 0.06 dB in all five cells, by two oracles — a streaming
implementation and a **grid-free** closed form (`M(t) = Σ w_j(t)·e^{jΩ(o_j − t)}` enumerated over the
exact 147 rational phases; no signal, no FFT, no bucketing). A time-stretched kernel misses that row by
1.2 dB at 22.5 kHz and 9.2 dB at 23.5 kHz, and would have given a round-trip delay of exactly 64.0
samples rather than 61.4.

### 6.1 Both axes — the modulation is gone, not merely smaller

Coherent carrier / worst phase / best phase, in dB, one round trip, block 64. The "was" column is the
cubic **re-measured by this same instrument**, not quoted:

| f, Hz | was: carrier / worst | now: carrier | now: worst | now: best |
|---|---|---|---|---|
| 5 000 | −0.05 / −0.08 | −0.000128 | −0.000392 | +0.000034 |
| 10 000 | −0.64 / −1.16 | −0.000054 | −0.000175 | +0.000025 |
| 15 000 | −2.59 / −5.14 | +0.000003 | −0.000053 | +0.000113 |
| **17 640** | **−4.17 / −9.27** | **+0.000160** | **+0.000016** | **+0.000385** |
| 19 000 | −4.98 / −12.20 | +0.000261 | +0.000151 | +0.000445 |
| 20 000 | −5.48 / −14.79 | −0.013301 | −0.013469 | −0.013113 |

Carrier, worst and best have **collapsed onto each other**: the modulation depth is at most 0.0004 dB
where it was up to 8.7. The residual is not noise and must not be asserted away — a Kaiser β = 8.6
window has a real passband ripple δ = 10^(−86.7/20) = 4.6e−5, i.e. ±0.0004 dB per stage and ±0.0008 dB
round trip, and the measured best phase (+0.00039 dB at 17.64 kHz) sits inside it. The suite's old
assertion `bestDb ≤ 0.0001` would have **failed the correct kernel**; it is now the derived bound.

The 20 kHz row is the designed band edge doing its job: 20 kHz is 0.907 of the 22.05 kHz Nyquist and
the cutoff is at 0.99, so −0.013 dB is the shoulder, not an error.

**And the shoulder past 20 kHz, which is the audible half of "brighter than the cubic" and was missing
from the first version of this table** (a crew round asked for it by name). Round trip, coherent
carrier:

| f | 20 000 | 20 500 | 21 000 | 21 500 | 22 000 |
|---|---|---|---|---|---|
| sinc | −0.013 | −0.327 | −1.962 | −6.770 | −20.43 |
| cubic | −5.475 | −5.667 | −5.813 | −5.906 | −5.941 |

(Both rows measured by the same instrument; the cubic's was not estimated.) Read it honestly: **the two
cross at about 21.2 kHz, and above that the sinc is the DARKER of the two**, because it
is doing the band-limiting the cubic simply refused to do. What the cubic passed up there was not signal
— it was the material that folded back and that a driven stage then demodulated (§4, §6.3). Everything
below 20 kHz is where "brighter" lives, and there the sinc is flat to 0.0002 dB against the cubic's
−4.17 at 17.64 kHz.

### 6.2 The decimating direction now has a stopband

48 → 44.1, tone above the 22.05 kHz output Nyquist, dB rms / sample peak:

| tone at 48 kHz | cubic | 64-tap sinc |
|---|---|---|
| 22 100 Hz | −3.0 / **0.0** | −9.08 / −9.08 |
| 22 500 Hz | −3.0 / **0.0** | −15.41 / −15.42 |
| 23 000 Hz | −3.1 / **0.0** | −27.45 / −27.45 |
| 23 500 Hz | −3.1 / **0.0** | −47.56 / −47.52 |
| 23 900 Hz | −3.1 / **0.0** | −88.77 / −85.92 |

**The 0 dB sample peak is gone**, and that is a separate fact from the rms: it existed because at phase
`t = 0` the cubic's weights were `(0,1,0,0)`, a bare sample pick with `|M(0)| = 1` at every frequency.
No phase of this kernel is a sample pick, so the peak now tracks the rms.

Where a 23 kHz tone lands, which used to be **two** strong components 1.7 dB apart covering
18.15–22.05 kHz:

| | cubic | sinc |
|---|---|---|
| 21 100 Hz (the `44100 − g` fold) | −5.33 | **−27.45** |
| 19 100 Hz (the folded image at `48000 − g`) | −7.04 | **−91.58** |

The second component — the one a crew round found and which is why the damaged band was not just the
top slice — is 84 dB below where it was.

### 6.3 The demodulation, which is the product half of this change

Measured through the real `NamStage` on real captures, 20 kHz tone at −18 dBFS, 44.1 kHz host, 2 s
coherent window. Both columns come from the SAME probe; the "before" column runs it against `main`
(52582a8) and reproduces §4's published figures exactly, which is what licenses the "after" column:

| `Luma Recti` (high gain) | before (cubic) | after (sinc) |
|---|---|---|
| carrier at 20 kHz | −32.16 dBFS | −35.81 dBFS |
| **the 100 Hz line it demodulates into** | **−17.65 dBFS** | **−96.54 dBFS** |
| …relative to its own carrier | **+14.51 dB** | **−60.73 dB** |

| `Angl Blackmore Clean` | before | after |
|---|---|---|
| carrier at 20 kHz | −13.41 dBFS | −9.56 dBFS |
| the 100 Hz line | −39.07 dBFS | **−117.36 dBFS** |
| …relative to its own carrier | −25.66 dB | −107.80 dB |

**The line dropped 78.9 dB on the high-gain capture and 78.3 dB on the clean one**, and the relation
that made P32 call this the loudest number in the task — a bass line *louder than the ultrasonic
carrier that produced it* — is inverted by 75 dB.

⚠️ **One honest qualification, because the band figure does NOT move by 79 dB.** The whole 20 Hz–4 kHz
energy on the 100 Hz lattice reads −17.20 dBFS before and −23.01 after: 5.8 dB, not 79. Broken down,
that band is **one bin** — 4 000 Hz at −23.0 dB — and the same bin reads −25.9 dB at a **48 kHz host,
where there is no resampler in the path at all**. It is the model's own 5th harmonic of 20 kHz folding
(5·20000 − 2·48000 = 4000 Hz), not anything the rate-match does. So the correct reading is: before, the
resampler put a term 5.8 dB **above** the model's own floor into the bottom of the spectrum; after, its
contribution has gone **under** that floor and what remains is the model being itself.

### 6.4 CPU — the price, on two machines

One mono round trip (both stages), 44.1 kHz host / 48 kHz model, best of five runs:

| machine | whole `NamStage` (model included) | cubic | 64-tap sinc | what the swap costs |
|---|---|---|---|---|
| arm64, Apple clang | 5.3–6.1 %RT mono | 0.0129 | **0.1447** | **+2.3 %** of the stage |
| x86-64 i9, gcc 14.2 | 11.0–11.6 %RT mono | 0.0525 | **0.3239** | **+2.4 %** of the stage |

Flat across host blocks 32 / 64 / 128 / 512 on both (spread under 2 %). Absolute %RT does not travel
between machines; **the fraction does**, and the two toolchains agree to 0.1 points.

🔴 **The inner loop's shape was chosen by measurement, and the obvious algebra is the slow one.**
Blending the two phase rows per TAP and summing once is identical arithmetic to two dot products
blended by a scalar, and looks more expensive (64 MACs + 64 lerps against 128 MACs). Measured at 64
taps on arm64: per-tap lerp **12.77 ns** per output (10.0 GMAC/s), two accumulators 16.70 ns, a
hand-written four-way split 21.21 ns. A float reduction is not associative so the compiler may not
restructure it: one dependency chain with a fused lerp vectorises, two compete, and splitting by hand
only spends registers. Picking the "cheaper" form would have cost 43 % more CPU.

### 6.5 Latency — derived from the geometry, then measured back

`reset()` primes `buf` with `kTaps` leading zeros and `pos = kHalf`, so buf[q] holds input sample
q − kTaps and output k is centred on input position **k·inPerOut − kHalf**. Every stage therefore
delays by exactly **kHalf = 32 of its own input samples** — the group delay of a symmetric 64-tap FIR.
The round trip is 32 host samples (down) plus 32 model samples (up) converted to host rate:

    D·(1 + hostSR/modelRunSR) = 32·(1 + 44100/48000) = 61.4000 host samples = 1.392 ms

**Measured back** from the carrier phase at 100 Hz and 200 Hz — below the first whole-period wrap, so
unambiguous — it reads **61.4000**, the geometry to four decimals. (At 500 Hz the same measurement
reads −26.8, which is 61.4 − 88.2 and equally true; that is the wrap the impulse onset resolves.)
Unlike the cubic, this kernel is symmetric and its phase delay is **flat with frequency**, so the old
"+0.018 samples at 10 kHz, +0.620 at 20 kHz" group-delay caveat died with the cubic.

Reported integers: **61** at 44.1 kHz, **96** at 96, **91** at 88.2, **53** at 32. `NamStage` no longer
restates the geometry at all — it asks `StreamResampler::delayInputSamples()`, because restating it is
exactly how the previous formula stayed 2.16 samples wrong through a release cycle.

⚠️ The exact-half rate MOVED with the kernel: with D = 2 the halves sat at h = 24000k − 36000 (12000,
36000, 60000 …); with D = 32 they sit at h = 750·(2k+1), and 60 kHz now gives a whole 72.0. A rounding
test that kept its old rate list would have gone silently blind.

### 6.6 Two contract changes that are NOT tuning

* **Identity ratio.** `sinc(0.99·n)` is not zero at integer n, so running an exactly-equal in/out rate
  through the general path would low-pass a caller who asked for no rate change. The class
  short-circuits exact identity to a **bit-exact copy** delayed by kHalf. The delay stays one formula
  for every ratio; only the filtering is skipped. `NamStage` never reaches this (its gate is 0.5 Hz).
* **The kernel is APPROXIMATING, not interpolating.** The cubic passed through its input samples
  (h(0) = 1, h(n≠0) = 0) and reproduced constant / linear / quadratic signals bit-exactly. This one
  does not: DC now lands within a derived 4.6e−6 (measured 3.6e−7, 46 % of samples still bit-exact),
  and a ramp within 4e−7 relative. Integer ratios are no longer sample-picking either — 96 → 48 kHz is
  a filtered decimation now, which is the point.

### 6.7 🔴 WHERE THE TRANSPARENCY CLAIM STOPS — it is a claim about the RATIO, not the kernel

Every other number in this document is measured at a 44.1 kHz host. The window is a fixed `kTaps = 64`
**input** samples, so its transition width in Hz scales with the **input** rate: at a 192 kHz host the
down leg is a 4:1 decimation and those same 64 taps buy a ~17 kHz transition, which starts inside the
audio band. Round trip host ↔ 48 kHz, coherent carrier in dB:

| host | 15 kHz | 17.64 kHz | 19 kHz | 20 kHz |
|---|---|---|---|---|
| 44 100 | +0.0000 | +0.0002 | +0.0003 | −0.0133 |
| 88 200 | −0.0003 | −0.0002 | −0.0005 | −0.0009 |
| 96 000 | +0.0003 | +0.0002 | −0.0004 | −0.0075 |
| **176 400** | +0.0001 | **−0.0374** | **−0.2352** | **−0.6147** |
| **192 000** | +0.0003 | **−0.0814** | **−0.3512** | **−0.7908** |

**So "transparent" is asserted for hosts up to 96 kHz and no further.** At 176.4 and 192 kHz this kernel
costs up to 0.8 dB in the top octave — still far better than the cubic, which had no anti-aliasing at
all and folded everything above 24 kHz back down, but not the flat response the rest of this document
describes. The rows are pinned in `felitronics_core_streamresampler_lptv_tests` in BOTH directions, so a
future kernel that fixed it has to come and edit them rather than quietly pass.

**The fix, if it is ever wanted, is a design change and not a tuning knob:** scale `kTaps` with
`max(1, inRate/outRate)`, which restores the transition width in Hz and costs CPU in exact proportion —
4× the taps at 192 kHz. That is a plan item, not a P34 one.

## 7. What is gated and what is only recorded here

`felitronics_core_streamresampler_lptv_tests` pins, so a kernel change has to come and edit this
document: the whole carrier / worst-phase / best-phase table for BOTH kernels, the period-147 closure,
block independence, the exact produced count, `produceExact`'s silence padding, all five decimation
rows, where a 23 kHz tone lands (both components), and the criterion — the round trip now adds
**−100.72 dBc** at 17.5 kHz where the cubic added −8.84.

🔴 **That suite carries a verbatim copy of the Catmull-Rom kernel, and it is load-bearing.** P32's
instrument proved itself by reading 0.00 dB on the unity ratio — an instrument that cannot read zero
cannot be trusted to read −4. That argument **dies** once the kernel under test is itself flat: a wrong
cutoff, a wrong tap centre, a mirrored table and an oracle bucketing by the wrong period all report a
beautiful 0.00 dB on a flat kernel. So liveness now runs both ways — the same instrument must first
reproduce the cubic's −4.17 / −9.27 before it is allowed to certify the sinc's zero. The minimality
half of the period-147 proof lives there for the same reason: on a flat kernel every divisor of 147 is
a period to within the noise, so the assertion would pass vacuously.

`felitronics_nam_tests` pins the same carriers through the real `NamStage` plumbing, the round-trip
delay measured two ways (impulse onset for the integer, carrier phase for the fraction), and the
48 kHz null — that a stage which has RUN the rate-matcher and been re-prepared at the model's own rate
renders bit-identically to one that never did.

Not gated, because it needs material that cannot ship in a repository — real captures and a real
performance: everything in §4 above and the demodulation table in §6.3. Those are reproducible from
the protocol but nothing will tell you when they rot; treat them as dated measurements, not invariants.

## 8. What `latencySamples()` is used for downstream — it is not only PDC

`nam::NamStage::latencySamples()` moved again with this kernel (4 → **61** at 44.1 kHz, 6 → **96** at
96, 6 → **91** at 88.2, 3 → **53** at 32), and inside this repository that number reaches nothing but a
report: `RigPlayer::latencySamples()` republishes the max of its two slots outward, and slot alignment
runs on `AlignmentTable::delayOf()` → `blendDelay()` / `lagTail_`, which never read it.

**Outside this repository it aligns audio.** OrbitCab delays its dry / bypass path by exactly this
number (`src/poweramp/PowerAmpRouter.cpp`, `src/core/CabEngine.cpp`), and orbit-amp does the same at
the dry end of its crossfade. P32 corrected the reported number from a guessed 6 to the true 3.84 and
moved the first comb notch of an on↔off crossfade from ~10.2 kHz out to ~136 kHz; the same arithmetic
holds here, with a residual of 0.40 samples at 44.1 kHz (notch at ~55 kHz) and 0.00 at 96.

🔴 **But the ABSOLUTE number is now 16× larger, and that is a product-visible change, not a rounding
one.** A host that reports 61 samples of PDC at 44.1 kHz instead of 4 will latency-compensate by
1.39 ms; anything that delays a dry path by this figure must be able to. Both shipped hosts also **sum
two rate-matching stages** for PDC (preamp + poweramp), so the number the host sees is **122 samples at
44.1 kHz** and 192 at 96.

### 8.1 What this moves in `orbit-amp`, measured rather than assumed

`orbit-amp` pins core **hybridly** (`CMakeLists.txt:69-80`): a sibling `../felitronics-core` checkout
wins over the release tag, so it meets a change like this **locally, at merge**, while its CI stays
green until the pin moves. Its whole test suite was built against this branch and run — seven targets,
zero compile errors, zero failures. **There is no orbit-amp test that breaks.**

What does break is not a test but a **silent clamp**:

> `orbit-amp/src/core/BypassWire.h:37` — `static constexpr int maxDelay = 64;`, justified at `:35-36`
> by *"The most a rate-match can cost: `ceil (3 * hostSR / modelSR) + 3`, which at 192 kHz against a
> 48 kHz pack is fifteen"*. The clamp bites at `:52` (`juce::jlimit (0, maxDelay, delay)`).

That formula is two generations stale — it is the guess P32 replaced with `2 + 2·h/m`, which P34
replaced with `D·(1 + h/m)`:

| host | as commented | v0.26.0 | **this branch** | what the wire delivers | short by |
|---|---|---|---|---|---|
| 44 100 | 6 | 4 | **61** | 61 | 0 |
| 48 001 | 7 | 4 | **64** | 64 | 0 |
| 88 200 | 9 | 6 | **91** | 64 | **27** |
| 96 000 | 9 | 6 | **96** | 64 | **32** |
| 192 000 | 15 | 10 | **160** | 64 | **96** |

So 44.1 kHz still fits, by three samples, and **every host above 48 kHz under-delays the bypass path
without complaining**. The consequence is exactly what that class exists to prevent — its own comment
calls blending against an undelayed copy "a comb" — and at 96 kHz the residual 32 samples put the first
null at 1500 Hz. The fix belongs in that repository and must be the GEOMETRY,
`D·(1 + hostSR/modelRunSR)` with `D` from `StreamResampler::delayInputSamples()`, not a bigger constant.

Two more addresses the merge touches, for whoever picks it up: `src/PluginProcessor.cpp:404` sums the
two blocks' latencies into the host's PDC (**8 → 122** samples at 44.1 kHz), and `:765` is where that
number is handed to the clamping wire (`:771`, `:785`, `:791`).

**OrbitCab is deliberately not listed.** It leaves with NAM and poweramp, and an earlier version of this
document pointed at three of its tests on the strength of P32's notes rather than a measurement.

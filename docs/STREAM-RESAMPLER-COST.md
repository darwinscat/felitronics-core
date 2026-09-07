# What `core::StreamResampler` costs — measured

`felitronics::core::StreamResampler` is a Catmull-Rom cubic interpolator with no anti-aliasing
filter of any kind. It rate-matches a live stream in and out of a rate-locked neural model
(`nam::NamStage` holds **two** instances per channel — `down` host→model before the network and
`up` model→host after it), and it engages only when `|hostSR − modelRunSR| > 0.5`
(`NamStage.cpp`). **At a 48 kHz host with a 48 kHz capture there is no resampler in the path at
all**, so everything below belongs to 44.1 kHz and other off-native rates.

The header used to justify the kernel with one sentence — *"the driven nonlinear stage masks the
interpolation images"* — and no number. This document is the number.

## 1. The mechanism: one thing, not two

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

## 2. Both axes, round trip 44.1 ↔ 48 kHz

Coherent carrier / worst phase, dB. Three independent oracles agree to 0.01 dB: an exact
rational-phase enumeration of the closed form, a brute-force simulation, and the shipped class driven
through `NamStage`'s own plumbing (cos and sin runs combined into a per-sample complex gain — the
class is linear, so no bucketing is needed and none is done).

| f, Hz | coherent | worst phase | best phase |
|---|---|---|---|
| 1 000 | −0.00 | −0.00 | −0.00 |
| 5 000 | −0.05 | −0.08 | −0.01 |
| 8 000 | −0.28 | −0.50 | −0.05 |
| 10 000 | −0.64 | −1.16 | −0.11 |
| 12 000 | −1.22 | −2.28 | −0.21 |
| 15 000 | −2.59 | −5.14 | −0.41 |
| **17 640** | **−4.17** | **−9.27** | −0.61 |
| 19 000 | −4.98 | −12.20 | −0.69 |
| 20 000 | −5.48 | −14.79 | −0.74 |

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

So the **worst-phase column is effectively a ceiling** (the shipped priming lands 0.02–0.9 dB off the
alignment-wide worst), while the **coherent carrier is an interference term** — stage 2 converts part of
stage 1's sidebands back onto the carrier, and how much depends on the alignment. The mean over alignments
is exactly the product of the two single-stage carriers (−5.06 dB at 17.64 kHz, −7.77 at 20 kHz), which is
what a "two independent stages" estimate would give and why that estimate is not the cascade.

## 3. The decimating direction has no stopband

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

dB rms / sample peak, unit input. The peak is an ENVELOPE figure — no single spectral line exceeds
−4.67 dB; it reaches 0 dB where the interpolation phases line up. What survives does not land in one place: a tone at *g* ∈ (22.05, 24) kHz
comes back as **two** strong components — `44100 − g` at about −5 dB and `g − 3900` at about −7 dB —
plus weaker terms near −45 dB, so the fold covers **18.15–22.05 kHz**, not just the top slice. Measured
on the shipped kernel, 2-second coherent window:

| tone in | components out |
|---|---|
| 22 100 Hz | 22 000 Hz −4.67 · 18 200 Hz −7.92 · 18 100 Hz −47.8 · 14 300 Hz −42.6 |
| 23 000 Hz | 21 100 Hz −5.33 · 19 100 Hz −7.04 · 17 200 Hz −45.9 · 15 200 Hz −43.2 |
| 23 900 Hz | 20 200 Hz −6.06 · 20 000 Hz −6.23 · 16 300 Hz −44.5 · 16 100 Hz −44.2 |

## 4. What the driven nonlinear stage actually does to it

Measured on the reference NAM DI through real captures, against an exact band-limited reference
conversion (`shipped − ideal`, both at 44.1 kHz):

* the **output** leg — the half with no possible causal masker, since nothing after
  `up.produceExact()` in `NamStage` or `RigPlayer` is nonlinear — is the **larger** half, by 5…10 dB
  in every run;
* against the **model's own aliasing floor** (what a WaveNet at 48 kHz folds down by itself, measured
  on an arm with no resampler in it), the **output leg alone** sits below it on a high-gain capture
  (3–24 dB) and **up to +9.8 dB above** it on a clean one from 15 kHz up. Counting the whole rate-match,
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
  fall on a 100 Hz grid, which is where they land. The clean capture does it too, 33 dB quieter
  (−39.1 dBFS at 100 Hz). A steady 20 kHz tone is not guitar, so read this as the mechanism rather than
  as a level; the level on real DI is the table above.

So the old one-line defence fails three ways: it never covered the carrier droop at all, it is
conditional where it does apply (below the model's own floor on a driven capture, up to +9.8 dB above
it on a clean one, output leg alone), and the driven stage does not mask the input leg's artifacts —
it moves them into the part of the spectrum where nothing masks anything.

## 5. If the kernel is ever changed

A polyphase windowed-sinc (Kaiser β = 8.6, `L` taps, `P` phases, linear interpolation between
phases) fits the same seam and the same RT contract (table built in `reset()`). Round trip, coherent
/ worst phase. At 32 and 64 taps the two columns **coincide** — the modulation is gone and only a
designed band edge remains; at 16 taps they still differ by 0.04–0.28 dB, i.e. it shortens the
modulation without removing it:

| kernel | 10 k | 15 k | 17.64 k | 19 k | 20 k | %RT (mono, both stages) | round-trip delay |
|---|---|---|---|---|---|---|---|
| shipped Catmull-Rom | −0.64/−1.16 | −2.59/−5.14 | −4.17/−9.27 | −4.98/−12.20 | −5.48/−14.79 | 0.0126 | **3.84 samples** |
| L=16 P=256 | −0.00/−0.00 | −0.20/−0.20 | −2.30/−2.34 | −5.30/−5.44 | −8.86/−9.14 | 0.0372 | 15.36 |
| L=32 P=256 | −0.00/−0.00 | 0.00/0.00 | **−0.11/−0.11** | −1.77/−1.77 | −6.10/−6.10 | 0.0790 | 30.71 |
| L=64 P=512 | 0.00/0.00 | 0.00/0.00 | 0.00/0.00 | −0.03/−0.03 | −2.60/−2.60 | 0.1451 | 61.41 |

For scale, the whole `NamStage` (model included) measures **5.3–6.1 %RT mono** at 44.1 kHz on an
Apple-Silicon build with a standard WaveNet capture, so a 32-tap pair costs about **+1.2 %** of what
the stage already spends. The price is latency: 30.7 host samples (0.70 ms) against today's 3.84
(0.087 ms). The 19–20 kHz roll-off of the sinc candidates is a chosen cutoff (0.94 of the lower
Nyquist), not a property of the kernel.

`rigplayer` consequences, checked in the code rather than assumed: slot alignment runs on
`AlignmentTable::delayOf()` → `blendDelay()` / `lagTail_`, and **none of them reads
`latencySamples()`** — `RigPlayer::latencySamples()` only republishes the max of the two slots
outward. So a kernel change moves the absolute PDC and nothing else inside the player, as long as
both slots take the same delay. `AlignmentTable::measureAlignment()` takes the rate as an argument
and defaults to 48 kHz, where the resampler is bypassed entirely; called at another rate it still
measures the same processed samples.

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

**These are ceilings, not sampled maxima.** At 44100/48000 = 147/160 the composite gain is periodic
with **exactly 147 output samples**: stage 2's position advances by 160/147 per output, so after 147
outputs it advances by exactly 160 — which is a whole number of stage-1 phase periods, since stage 1
advances by 147/160 and repeats every 160. Every phase pair that the topology can produce therefore
occurs inside one period of 147, and min/max over it are exact. (Checked numerically: the deviation
from period-147 repetition is 2.6e-13 over five periods, no proper divisor of 147 is a period, and
the statistics over 23 520 samples are identical.)

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

dB rms / peak, unit input. Everything that survives lands in 20.1–22.05 kHz.

## 4. What the driven nonlinear stage actually does to it

Measured on the reference NAM DI through real captures, against an exact band-limited reference
conversion (`shipped − ideal`, both at 44.1 kHz):

* the **output** leg — the half with no possible causal masker, since nothing after
  `up.produceExact()` in `NamStage` or `RigPlayer` is nonlinear — is the **larger** half, by 5…10 dB
  in every run;
* against the **model's own aliasing floor** (what a WaveNet at 48 kHz folds down by itself, measured
  on an arm with no resampler in it), the added artifacts sit **below** it on a high-gain capture and
  **up to +9.8 dB above** it on a clean one, from 15 kHz up;
* driving harder makes it **worse**: over a 42 dB input sweep the total error rose +11.4 dB relative
  to the program, the folded-alias term +12.5 dB, and the model's own >22.05 kHz output +13.2 dB —
  because "driven" is precisely what fills the band the un-filtered decimation folds back.

So the old one-line defence is conditional (true for a heavily driven capture, false for a clean
one), it never covered the carrier droop at all, and the variable it invokes moves the wrong way.

## 5. If the kernel is ever changed

A polyphase windowed-sinc (Kaiser β = 8.6, `L` taps, `P` phases, linear interpolation between
phases) fits the same seam and the same RT contract (table built in `reset()`). Round trip, coherent
/ worst phase — note that for these the two columns **coincide**, i.e. the modulation is gone and
only a designed band edge remains:

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

`rigplayer` consequences, checked: both slots would take the same new delay, so their **relative**
alignment does not move — only the absolute PDC; `AlignmentTable` measures at 48 kHz by default,
where the resampler is bypassed entirely; `nam::blendDelay` / `lagTail_` is an integer line unaffected
by a delay common to both slots.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### oversampling · saturation · limiter · poweramp — a strict oversampler that is flat to 20 kHz, as an option

`PolyphaseOversampler` centres a fixed cutoff at 0.45 fs, and at 44.1 kHz that sits below the top of the
audio band: **one round trip costs −1.80 dB at 19 kHz and −15.55 dB at 20 kHz**, and more taps make
20 kHz worse (−13.71 at 32, −19.17 at 120). The cutoff is not an oversight, and its derivation is now in
the header: the class wraps a nonlinearity, so an image of content between 20 kHz and fs/2 that the
interpolator only partly rejects is multiplied with that content and lands in the audio band. The
design is therefore **strict** — the transition must finish below fs/2.

- **`oversampling::CascadeOversampler`** keeps the guard and moves the band edge: a Kaiser 2x stage whose
  length and cutoff follow from the sample rate by a measured rule, then halfband 2x stages. Images and
  aliases of anything below fs/2 at **−91 dB or lower**, one pass within **0.0043 dB** up to 20 kHz, over
  14 rates from 8 kHz to 768 kHz and factors 2–16. Powers of two only; refuses (law 11b) what it cannot
  build.
- **The price is latency, not CPU**: **131 base samples at 44.1 kHz 4x** (63 for the Kaiser stage) for
  about the same multiply count; **76 at 48 kHz, 28 at 88.2 kHz**, where it is also cheaper than the
  Kaiser stage.
- **A halfband first stage — the cheaper-looking cascade — cannot do this**, at any length:
  H(f) + H(Fs/2 − f) = 1 makes its image rejection at fs − a equal its pass-band deviation at a. A 79 + 23
  pair that loses 0.41 dB at 20 kHz leaves that tone's image at −32.5 dB; one flat to 20 kHz still passes
  don't-care content's images at −58 / −35 / −22 dB (r = 0.46 / 0.47 / 0.48) — both pinned in
  `felitronics_oversampling_tests`. Through a single-ended tube stage at +12 dB that is −48 … −56 dBFS of
  intermodulation in the audio band (its 2nd-order line at 0.9–1.8 kHz) from −20 dBFS of content at
  0.48–0.49 fs, against under −150 for the strict designs (P31 stand, not a test).
- **`oversampling::Topology`** (`Kaiser` | `Cascade`) and `oversampling::Oversampler`, which is either.
  `Saturator::prepare`, `TruePeakLimiterConfig::topology` and `PowerAmpStage::prepare` take it; **the
  default is `Kaiser` everywhere, and under it nothing changes** — same refusals, same latency, same bits
  (the whole suite passes unchanged). Under `Cascade`: `tapsPerPhase` is still range-checked and
  otherwise unused; the Saturator's dry path is delayed by the cascade's round trip; `PowerAmpStage`
  rounds the factor down to a power of two and clamps the design rate, as it clamps everything.
- **The limiter's ceiling under the cascade**, measured with the ceiling suite's own witnesses: at
  44.1 kHz 4fs/9 is now delivered flat, so the grid allowance at 8x rises from 0.108 to **0.133 dB**
  (16x: 0.027 → 0.033); the 1.15 dB modulation envelope holds (worst 1.344 dB at 2x, inside 2.399), and at
  4x and 8x the dense excess is lower than Kaiser's (0.49 / 0.48 against 0.81 / 0.75 dB).
- **Unchanged on purpose**: `ReferenceTruePeakMeter` stays `PolyphaseOversampler` at 4x/32 — it is the
  certified unit — and no stage's default moved. `MasteringChain` and the wasm ABI do not offer the
  topology.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### convolution — 🔴 SOUND CHANGE: an un-normalized IR keeps its level across sample rates

`CabConvolver` prepared with `normalize=false` convolved **louder in proportion to the host's sample
rate**. A 48 kHz IR played **+6.02 dB** hotter in a 96 kHz session than in a 48 kHz one, **+5.28** at
88.2, **+12.04** at 192, and **−0.74** at 44.1 — the same file, the same knob, a different clock. It now
plays at the same level at every rate, which **moves shipped sound** by exactly those figures. **Nothing
changes for a load that is not resampled** — the IR already at the host's rate (within the 1 ppm match
tolerance), or one whose rate is unknown — where the multiplier is exactly 1 and the staged taps are
byte-identical to what v0.33.0 published. A 48 kHz host is NOT by itself such a case: a 96 kHz IR played
6.02 dB quiet there and now does not.

- **Why it happened.** `resampleIr` normalizes each output tap to unity DC, which preserves a WAVEFORM's
  amplitude — right for a signal, and what the function is for. A convolution's gain is not an amplitude:
  it is a sum over taps, so it is proportional to how many of them fit into a second. Resampling 48 → 96 kHz
  produces twice as many taps of the same amplitude, hence twice the gain. JUCE's convolver applies exactly
  this factor on exactly this branch — `juce_Convolution.cpp:790-792` in 8.0.14: `normalise == yes` measures
  the resampled taps, `else` calls `resampled.applyGain (originalSampleRate / processSpec.sampleRate)`. The
  interim JUCE backend this replaced therefore had it; the replacement carried over the normalized half of
  that `if` and not the other one.
- **And this repo's own spec said so.** `docs/migration/convolution-core-api.md:30` — the parity table whose
  file header calls it "the authoritative spec … without changing the sound" — has carried the row
  **`Normalise::no still applies a gain | makeEngine, applyGain(irSr/hostSr) (:792) | after any resample,
  multiply IR by irSr/hostSr`** since the migration was written. The requirement was not missed for want of
  being written down: nothing ever checked a row of that table, and the same file's pseudocode ten lines
  later contradicted it ("every KNOWN rate" where the table says "after any resample"). The pseudocode is
  corrected here; the table was right all along. And the golden-test tier that spec defines for a
  non-host-rate IR is `Normalise::yes` — the one tier that structurally cannot show this.
- **Why only this path.** With `normalize=true` the reference-unity gain is measured from the FINAL,
  resampled taps, so it absorbs the factor whole. The defect was visible only on the one path nobody
  measured that way.
- **The fix** is one factor, `irSr/hostSr`, applied where the loader already applies its one gain. It is
  named — `convolution::convolutionRateGain(inSr, outSr)`, beside the resampler — because a second copy of
  this arithmetic was already written by hand in a shipped product, which is the class of defect where two
  copies must agree forever. A caller that resamples an IR itself should use it instead of spelling it again.
  Rates it cannot make a factor from (a NaN, a zero, a negative, an infinity, or a finite pair whose quotient
  overflows) return 1 — the same rates `resampleIr` refuses get no compensation.
- **The factor is NOT clamped, and its sibling is.** The normalization gain is measured from the IR's own
  content, so ±30 dB stops a pathological IR blasting; the rate factor is arithmetic on two numbers the
  caller supplied, and clamping it would quietly deliver a different filter than the caller asked for. It
  spans 4.7e−10 to 4.3e9 (−186 to **+192.6 dB**), which is what `resampleIr`'s length and position gates
  leave reachable. A WAV's rate is a `uint32` and nothing in this family validates it, so a garbage-but-
  **finite** header is the one broken metadata a real file can carry — and it now plays LOUD where it used
  to play quiet: a file claiming 352800 Hz on a 48 kHz host is +17.3 dB, one claiming 5e6 Hz is +40.4 dB.
  Both are correct by this contract and both are garbage. A consumer that loads untrusted files should
  bound the rate before this, the way orbit-amp's loader already refuses anything outside 8 kHz…768 kHz.
- **What it does NOT fix.** The density term is all it removes. A resample still drops the kernel's
  pre-ringing that would fall before output sample 0, and what that costs depends on the onset, not on this
  factor — and it is a RIPPLE, not a loss: the kernel's nearest pre-ring lobes are negative, so cutting
  them ADDS level. Measured on a lone impulse at input index `lead`, 96 → 44.1 kHz: **−2.57 dB** at lead 0,
  **+0.65** at 1, **+0.95** at 2, −0.05 at 3, −0.61 at 4, and so on down to exactly nothing from `halfTaps`
  (32 input samples) of lead onward, which is the window's own backward reach. A smoother onset pays far
  less: a one-pole `exp(−n/τ)` at full scale with no lead loses 0.21 dB at 44.1 kHz for τ = 1 input sample
  and 0.004 dB for τ = 1 ms, which is what a real decay looks like. It is older than this change and
  untouched by it, but it is the one case where an **onset-trimmed** IR's level still moves with the rate —
  so a measurement that trims the front must not charge it here **in either direction**.
- **A diagnostic that had gone stale.** `irNormalizationGainDb()` floored its argument at 1e−6, which sat
  below everything the old code could apply (the normalization clamps at −30 dB) and above what this one
  can: an IR file claiming 0.024 Hz on a 48 kHz host applies 5.0e−7 and the reading said −120.0000 dB for a
  gain that is −126.0206. The floor is now below anything the loader can produce — **and a NaN gain now
  reports a NaN**, which is the one reading that did move on the normalized path: `std::max(a, b)` is
  `(a < b) ? b : a` and every comparison against a NaN is false, so `max(floor, NaN)` returned the FLOOR,
  a plausible −120.0000 dB for a gain that is not a number (a NaN tap in a float32 WAV is enough). The
  linear accessor was right throughout.
- **`irNormalizationGain()` now reports whichever gain was applied** — reference-unity on the normalized
  path, the rate factor on the verbatim one — and is still exactly `1.0f` when nothing was applied. The name
  is historical; the contract is "what was multiplied in", which is what a reference nulling against the
  engine needs.
- **One deliberate departure from the JUCE branch cited above:** JUCE applies the factor even when its
  resampler early-returned unchanged (equal rates), so it multiplies by a number that is only approximately
  1. This applies it only after a real resample, which is what keeps the byte-verbatim guarantee at matched
  rates. The two differ by at most `kRateMatchTolerance` relative — **8.7e−6 dB** — and the spec table's own
  wording ("after any resample") is the one this follows.
- **Pinned by** a property over eight rate pairs × five physical frequencies (44.1 / 48 / 88.2 / 96 / 192 kHz
  hosts, 44.1 / 48 / 96 kHz sources), read by a DTFT that shares no code with the loader, at a tolerance of
  0.01 dB — **13×** the 64-tap Kaiser kernel's own passband ripple (β = 8 → δ = 8.6e−5 → 0.00075 dB), **18×**
  the worst reading over the grid (5.65e−4 dB), and **74× below** the 0.736 dB the weakest cell moves without
  the fix. The kernel's column sums — the quantity the factor actually inverts — were measured across 64
  phases of every standard rate pair and stray at most 1.6e−4 dB from the ratio, which is the part of that
  argument that is a measurement rather than an engineering formula. Fifteen mutations of the change are each
  caught, among them "remove it", "invert it", "apply it only when upsampling", "apply it only to mono",
  "apply it on the normalized path too", "pre-scale the taps before the normalization is measured" and
  "delete the resample and keep the gain" — that last one is why an existing P67 assertion was rewritten:
  "every tap moved" stopped witnessing that the resampler ran the moment a gain could move every tap by
  itself, so the witness is now "no single scalar maps the input onto these taps".

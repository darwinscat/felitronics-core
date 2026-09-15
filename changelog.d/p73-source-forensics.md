### Added
- **`analysis::SourceForensics` — what a delivered file actually was, as far as the samples can prove it.**
  Two families of evidence and no attribution: a long-term (Welch) **spectral wall** — where the upper band
  limit is, in Hz *and* as a fraction of this rate's Nyquist, how deep the drop is, how wide the transition
  is, how far the spectrum comes back above it, and how far up it is empty (the upsampling tell) — and the
  **sample grid**: the coarsest dyadic grid `2^-k` every finite sample lies on, the shortest normalised PCM
  word that holds them all, the distribution of `k` behind that maximum, and the count of distinct sample
  values. Composed over the `analysis::SpectrumFrames` producer; `setParams` / `prepare` / `process` /
  `finish` / `reset`, `Storage` + `storageFor()` published before the allocation (law 11d), `maxBlock` sizes
  nothing, and every published field is bit-identical under arbitrary re-slicing of the stream into
  `process()` calls (law 8a).
- **`fcore_measure forensics`** prints the whole report, every float as a raw IEEE-754 bit pattern.

### Notes
- **A band limit is a suffix property, and its FLOOR is read past the transition.** The strict suffix
  maximum is anchored at the transition's end, not at the winning boundary: the winner sits on the cell
  that contains the band limit — that is what makes it the winner — so anchored there it reports a clean
  wall with 80 dB of "recovery" above it and a strict drop near zero (measured: 805 of 2350 constructed
  cutoffs, worst 80.4 dB, and one NEGATIVE drop).
- **Is there a floor to reach at all?** The floor reference is the median of the span above the edge, so on
  a stopband that is still descending it sits halfway down the descent and the transition "ends" half a
  span past the edge: an identical 16 kHz brickwall measured 49.8 Hz of transition over a flat floor and
  1594 Hz over one decaying at 15 dB/kHz, with `sharp` flipping false at 10 dB/kHz and true again at 25.
  The span's two halves are now compared; disagreement means there is no single floor, the width is a lower
  bound, and `transitionClipped` says so — a flag that was otherwise provably unreachable, so the promise
  attached to it had never been kept.
- **The grid's exact reading is a maximum, so a robust one is published beside it.** On float-rendered
  material — a 16-bit programme with a float fade-out, which is every real render — the fade's samples take
  `gridExponent` to 85 or 149 and the exact word length goes dark. `robustPcmBits` is the shortest word
  holding all but `gridOutlierFraction` (5 % by default, set by the length of a real fade) of the non-zero
  samples, computed from the histogram already accumulated; a tolerance of 0 makes it the exact reading.
- **A geometry that cannot finish is refused, not discovered at finish().** The edge search sorts one
  plateau span per candidate, so `fftOrder 22` with 0.01 Hz cells and a 12 kHz plateau span is 2.1 million
  candidates sorting a million doubles each. `storageFor()` and `prepare()` refuse it; the default geometry
  is 19240.
- **A band limit is a SUFFIX property.** The obvious construction — the steepest local descent, then the
  level just above it — reports a deep notch (a band-stop, a comb, a room null) as a 100 dB wall with a
  one-cell transition while full-power spectrum resumes 400 Hz higher, and it can select that notch and
  thereby discard a genuine wall above it. So the floor of the drop is the loudest cell *anywhere* above the
  candidate, and the search maximises that conservative drop. A monotone roll-off then fails on the drop
  itself (~10 dB) rather than on a transition-width technicality, and both negative cases are handled by one
  mechanism.
- **An edge that is not there is not given a frequency.** An argmax always returns something — the lowest
  candidate on a falling spectrum, the highest on a flat one — so `minDropDb` gates whether an edge exists
  at all: below it the report is `valid = false, ShallowerThanMinDrop`, with every number still published as
  evidence.
- **A narrow line above the edge is forgiven to a published rank and no further.** A 19 kHz whine or a
  leaked pilot tone is the loudest thing above a real wall; the search therefore runs on a 3-cell
  median-filtered copy of the cells and takes the `(t+1)`-th largest of the suffix, `t` capped by
  `exemptCells` and by a tenth of the suffix, so the exemption self-disables near Nyquist. What was forgiven
  is published (`sufMaxPower`, `strictDropDb`, `exemptedCells`, `recoveryDb`), and the emptiness test runs
  on per-bin maxima, which see every line.
- **The position near Nyquist is not attributed.** A converter's anti-alias filter, a 320 kbit/s codec and a
  genuinely band-limited master all live at 0.9-0.95 of Nyquist; `nearNyquist` marks that band and the
  instrument stops there. `nearNyquist == false` is not a claim that anything compressed the file.
- **Two edges, and which one is the report.** The primary is whichever candidate has the stronger
  conservative drop, so a 16 kHz codec wall inside a 20.5 kHz export filter reports the inner edge as the
  primary and the outer as the second — but only when the shelf between them is deep. With a shallow shelf
  the outer edge wins instead and the inner one is not reported at all: a second search BELOW the primary
  is structurally useless, since everything above such a candidate includes the primary's own plateau.
- **The PCM range is part of the claim.** A b-bit word carries `i/2^(b-1)` for `i` in
  `[-2^(b-1), 2^(b-1) - 1]`, so `-1.0` is a PCM sample and `+1.0` is not, at any depth. A sample outside
  `[-1, +1)` therefore withholds `minExactPcmBits` (reason `OutsidePcmRange`) while `gridExponent` stays
  valid, and the signed minimum and maximum the claim rests on are published.
- **What the grid proves, and in which direction.** `minExactPcmBits = k + 1` is the *shortest* normalised
  PCM word that holds the observed samples exactly; it does not bound the source's word length from above (a
  24-bit file carrying a padded 16-bit master is indistinguishable from a 16-bit one), and the *absence* of
  zero low bits proves nothing at all — 16-bit dithered up into 24 reports `<= 24`, never 16. A grid finer
  than `2^-23` is refused as a PCM word rather than reported as 32: float32's top binade is spaced `2^-24`,
  so a 32-bit stream reads `k = 24` near full scale and `k = 31` lower down, while a 32-bit *container*
  carrying 24-bit content reads 24 and one carrying 16-bit content reads 16.

### `convolution` — the IR resampler stops inventing a cabinet's top octave; a one-tap IR and a nearly-equal rate load

`convolution::resampleIr` treated what lies past either end of an impulse response as the weighted mean of the
samples it had, not as silence: a tap outside the input was skipped before its weight was added, so every edge
sample was divided by only the part of its window that landed on the input. A cabinet IR starts at its onset, and
`CabConvolver` resamples whenever the file's rate is not the host's — the ordinary case — so a 48 kHz cabinet in a
44.1 kHz session played a broadband floor over its own top octave. Twenty-one factory cabinets, measured outside the
tree, worst 1/6-octave band against each cabinet's own response: at 44.1 kHz +7.2 dB at 16 kHz and +25.1 dB at
18 kHz, now +0.08 and -0.45 dB; at 96 kHz +18.5 dB at 18 kHz, now -0.12 dB.

- **Zeros outside the input.** Every output sample divides by the weight of its whole window; a tap past either end
  adds weight, not signal. Resampling `[zeros, x, zeros]` now equals resampling `x` shifted by whole samples. What
  zeros cannot give back is the kernel's pre-ringing that would fall before output sample 0 — an IR that starts at
  sample 0 keeps it only by adding delay. It is a cut of at most 0.45 dB at 18 kHz on the factory set and grows with
  how abruptly an IR starts: the test's cabinet-like fixture reads +1.23 dB at 20 kHz at 96 kHz.
- **At least one sample.** The length is `inLen * ratio` rounded but never zero, as JUCE's
  `resampleImpulseResponse` had it. A one-tap IR at 96 -> 44.1 kHz rounded to nothing, and `CabConvolver` then
  published nothing without a word: the previous cabinet kept playing. The result is empty now only for a rate that
  is not a positive finite number, a length past `INT_MAX` — refused before the cast, where 2^32 + 1 used to wrap to
  ONE sample — or a ratio so small that output sample 0's position is past `INT_MAX` (the one-sample floor is what
  made that `(int) floor(t)` reachable; it is refused rather than undefined).
- **`CabConvolver::kRateMatchTolerance`** (new, relative `1e-6`, the same as orbit-amp's `sameRate`). Two rates this
  close are one rate and the IR loads verbatim; the exact comparison sent a host reporting 48000.0000001 through the
  resampler, which band-limits at 0.95 of Nyquist and moves every tap. Measured outside the tree on two factory
  cabinets, against a 2048-tap resample: played verbatim 1 ppm off its rate, one reads -87.5 dB over its first 100
  ms and the other -72 dB over its first second, where running the 64-tap resampler costs -73 and -81 dB; at 10 ppm
  verbatim is the worse of the two.
- **A load that stages nothing is ignored whole.** `buildAndStage` overwrote the retained taps before it knew the
  load was empty, and returned after. With a load still pending from mid-crossfade, the retry then used the old
  length on an emptied or narrowed store: a read past the end of a vector (a zero-length load, or a mono load over a
  pending stereo one — ASan container-overflow, libc++ hardening abort), or, with a null data pointer, a retry
  refused forever, `isBusy()` stuck true and neither the pending IR nor the new one ever reaching the convolver. The
  taps, their gain and the pending geometry are now staged in locals and committed together; a zero or negative
  length, a null channel array or plane, or a rate the resampler refuses leaves the playing IR, `stagedTaps()` and
  any pending retry exactly as they were.
  - **Behaviour:** after such a load `stagedTaps()` still holds the previous taps (it used to be emptied) — which are
    the ones playing. A null plane with a positive length is ignored instead of dereferenced.
- **Tests.** `felitronics_convolution_resampler_tests`: shift invariance at both edges over ten rate and radius cases
  and five inputs (in front to 1e-6, behind to the bit — on `main` 1437 and 1053 misses, 3.25e-2 worst on the
  cabinet fixture), the edge taps against an independent long-double recomputation of the specification, a
  hand-worked 105/104, a cabinet-like IR's bands against its own response and against the untruncated resample, the
  one-sample floor at seven short IRs, and the refusals. `felitronics_convolution_cabconvolver_tests`: the
  tolerance witness (verbatim to the bit at 48000.0000001 both ways and at half the tolerance at 192 kHz, resampled at
  twice it both ways), a one-tap IR off-rate staging, publishing and playing at four rate pairs, and six loads that
  stage nothing over a pending one plus a one-tap load that wins as the latest.

### deesser · dynamiceq — `detectorLevelDb()`: the level the detector sees, for drawing it against the threshold

`DynamicEqBand::detectorLevelDb()` and `DeEsser::detectorLevelDb()` return the linked sidechain envelope in dB: the
exact number the gain computer was last handed, so a UI can draw a live level line against the threshold, with
`gainReductionDb()` / `dynamicDeltaDb()` as the dip below it. Fed through the same curve and ballistics, the readout
reproduces the module's own gain bit for bit, sample by sample, audio and pauses alike.

- It is written where the curve's argument is made, so every step that advances the envelope writes it: the audio
  loop and a pause (`process` with zero channels) alike. `reset()`, `prepare()` and a mode or listen switch return
  it to the 1e-9 level floor.
- The floor is -180 dB exactly on `DynamicEqBand`, which narrows the dB to float before the curve. The de-esser's
  split-band detector does not narrow, so its floor is -180.00000025 dB.
- `DeEsser` reports the detector in the signal path. In `DynamicEq` mode with `listen` off that is the
  `DynamicEqBand` delegate's. In `SplitBand`, or with `listen` on, it is the local one.
- It is always finite: the sidechain input is gated (non-finite -> 0, |x| clamped to 1e6) before the band-pass. A
  grid of NaN, ±inf and ±FLT_MAX input reads at most 118.2 dB.
- It is a plain member read like `gainReductionDb()`: no allocation, no lock. Audio and meters are bit-identical
  to before.

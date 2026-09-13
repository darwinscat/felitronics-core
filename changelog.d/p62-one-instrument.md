### `analysis` · `mastering` · `tools` — one true-peak instrument aims and certifies a delivered ceiling

`TargetLoudnessSolver` promised a delivered file `<= maxTruePeakDbTp`, judged every render with
`analysis::TruePeakMeter`, and the file was certified by `fcore_measure`'s reference filter — two different
instruments. On bright material, and at the high delivery rates where the short filter stops interpolating, the
first read under the second by more than the solver's whole 0.05 dB aim, so a render the solver called feasible
was delivered over its promise. Fourteen real programmes delivered at six rates, outside the tree: 12 of 84
deliveries certified above -1 dBTP before, 0 after.

- **`analysis::ReferenceTruePeakMeter`** (new). The reference true peak as a module class: `PolyphaseOversampler`
  at 4x / 32 taps per phase at every rate, the maximum of |x| over the 4x stream floored at the sample peak, and
  `drain()` for the FIR's tail. Law 11 in the house order; a channel that stops is DRAINED at that moment rather than
  dropped — its pending peak was submitted and belongs to the reading, and it returns from silence — so a zero-width
  call is a pause that drains every channel. `prepare()` refuses a non-positive block as well as a bad rate or width;
  `storageFor (rate, maxBlock, channels)` equals what it allocates; nothing is allocated in `process()`/`drain()`.
  `analysis` now links `oversampling` (which depends only on `core`).
- **`fcore::Probe` measures through it.** `fcore_measure` and the browser's `fc_probe` report the same numbers,
  byte for byte (`--precise`, 40 runs over a real corpus against the previous binary). One sequence behaves
  differently, and neither shipped caller makes it: a channel left out of a narrower call and then given audio again
  no longer replays its pre-gap history. A narrowing stream that simply ends reads what it read before, to the bit.
- **The solver reads every render with the reference.** The ceiling is aimed, feasibility judged and
  `measured.truePeakDbTp` reported on the certificate's arithmetic: the reported number is now the certificate of
  the delivered samples, bit for bit (above the dB floors — silence is still spelled -200 dB, from the same float
  threshold as before).
  - **Behaviour:** a render that only the old meter called feasible now costs one more pass (10 extra passes over
    the 84 real deliveries), and the loudness it reaches is unchanged (worst move 0.0004 LU); a pass costs 10-15 %
    more (the reference filter on the whole programme, 60 s stereo, 44.1 -> 192 kHz).
  - **Budgets:** `solveBytes()` — and so `FC_NEED_SOLVE` — is the loudness meter plus 21 008 B for a stereo
    reference meter at any rate, where it was the loudness meter plus 392-296-248 B of `TruePeakMeter` plus a
    512 B drain buffer. `TargetLoudnessSolver::kDrainFrames` is gone: nothing is drained from a buffer any more.
- **`felitronics_truepeak_instrument_gap_tests`** (new) owns the number "how far apart the two meters read": five
  materials (music, drums, bright noise, a 16 kHz burst at its worst phase, a click at its worst offset) delivered
  at 44.1-192 kHz through the chain and pinned per cell. The worst is the drums at 44.1 kHz, where the cheap meter
  reads 0.2995 dB under the reference and 0.3267 dB under the band-limited truth (read at a continuous time, not on
  a zero-padding grid); the burst rows are the grid's closed form at 96 and 192 kHz, and at 176.4 kHz and above the
  cheap meter is shown to be a sample-peak meter. The reference is not the truth either, and that is pinned too:
  0.0272 dB under it on the drums, 0.3270 dB under it on a click flat to 0.45 fs at 48 kHz.
- **`felitronics_delivered_ceiling_tests`** (new): every source rate to every delivery rate of the six solves,
  certifies at or under the promise, and reports the certificate exactly.

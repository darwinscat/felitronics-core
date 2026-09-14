### analysis — `HumDetector`: mains hum, and the notes that look like it

`analysis::HumDetector` measures mains hum: a narrow, stationary line at 50 or 60 Hz with a comb of exact
multiples, looked for only inside the QUIET stretches of a programme, because music masks it. It reports the
line's interpolated position in Hz, its level over the local background in dB, which harmonics were found,
and which stretches were used — in sample coordinates. It runs on `analysis::SpectrumFrames`, is bit-identical
under arbitrary re-slicing of the stream into `process()` calls (law 8a), and publishes its whole budget
through `storageFor()` before allocating (law 11d). `fcore_measure hum` prints the report, every float as a
raw bit pattern.

**No answer reads as "clean" when the instrument could not look.** Seven named incompletenesses, not one of
them a zero that could be mistaken for absence: no quiet stretch at all; exactly one (because "stands still
between stretches" cannot be tested inside one); fewer than two stretches long enough to show the line twice;
a window too short to separate 49.0 Hz from 50.0 Hz; every frame holed; a mains-compatible line that was
measured and did NOT stand still; and a comb of harmonics whose base lies outside this scope. The evidence —
positions, levels, prominences, counts, the strongest peak of each search window whether accepted or not —
is published in every case.

What separates hum from a bass note (G1 is 1.0 Hz from the mains; A#1 and B1 are 1.7 Hz from 60) is the
position measured INSIDE the bin by a three-bin parabola, a 0.5 Hz tolerance that a drifting grid fits and a
note cannot, and the requirement that the same line stand still across at least two quiet stretches and
across every frame of them. The comb is reported, not required: a bass guitar's partials are near-exact
multiples too, so gating on a harmonic count would cost a false clean and buy little.

Four things in it are not obvious and are the reason it works:

* **Every local maximum of the search window is examined, not the strongest one.** A bass note 1 Hz away and
  14 dB louder owns the window's argmax; reading only that would lose the hum underneath it, with the
  resolution paid for and unused. A lone Hann-windowed tone's sampled skirt is monotone, so enumerating
  costs nothing in false positives.
* **Resolution is a duration, stated in bins.** Separating two Hann main lobes 1.0 Hz apart needs 2.7 bins
  between them — 2.0 bins is where the dip disappears — so the bin must be 0.37 Hz or finer: N = 2^17 at
  48 kHz, a 2.73 s window. `fftOrder = 0` picks the shortest window that achieves it at the file's rate.
* **A line must be prominent AND loud.** A ratio alone certifies arithmetic: a programme that is one pure
  tone has a spectrum of 1e-23 elsewhere, and a maximum of that residue stands 18–37 dB over the residue
  beside it. `minLevelDbfs` (−100 dBFS) is the other half of the test. A frame of exact digital silence is
  excluded for the same reason — it cannot have seen anything, so it is not evidence of absence.
* **A quiet stretch must show the line twice.** One frame is one periodogram, whose tail crosses a 10 dB gate
  often enough that a looped or duplicated quiet passage repeated the excursion and passed stationarity with
  it — measured, 75 of 1200 noise-only files and 68 of 1200 dither-only files. Two frames per stretch takes
  both to zero and costs no sensitivity.

The quiet gate measures the programme with the candidate bands removed, so a hum loud enough to fail the gate
cannot censor its own detection.

`mains` means a mains-COMPATIBLE stationary line was measured, not a causal claim: a synthesised 50.000 Hz
pedal with exact harmonics can be sample-for-sample what an interference pickup leaves. Four limits no
threshold can remove are named in the header, each with the field a consumer reads instead.

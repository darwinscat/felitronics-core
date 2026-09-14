### analysis — `LowEnd`: wide bass, and which note owns the bottom

New offline instrument `analysis::LowEnd`, for a master headed to a lacquer. A cutter head writes the
mono sum laterally and the difference vertically, so out-of-phase low end is physics rather than taste.
Two halves, one clock:

- **Wide bass.** An LR4 split at `crossoverHz` (default 120 Hz, `eq::Crossover2`), Mid/Side, and the
  **side energy fraction** `S/(M+S)` of the low band — over absolute 10 ms blocks (the same
  `lround(0.01·fs)` grid `analysis::LoudnessMeter` uses, so the two instruments name the same
  intervals), as a duration-weighted 100-bin distribution, and integrated. A fraction and not the
  requested `S/M` ratio because `L = −R` is a real master that makes Mid *exactly* zero; both raw
  energies are published, so `S/M = f/(1−f)` is one line away. Mono reads exactly 0 and is **valid**;
  digital silence is the 0/0 and reads `NoEnergy` rather than a zero that would look mono.
  Three extremum coordinates, because they are three different questions: the worst fraction (with its
  energy beside it, so an accidental 1.0 can be weighed), the loudest block, and the greatest **vertical
  modulation** — the one a cutting engineer asks for first and which neither of the others identifies.
  The high band's own Mid/Side pair and the unfiltered pair are published too; low + high is an LR4
  allpass, not the input, and the header says so instead of implying additivity.
- **The dominant low note.** 30–300 Hz folded to semitone bands (40 of them at A4 = 440), from the
  shared `analysis::SpectrumFrames`, with the repository's fractional-edge power integration — but
  computed from a precomputed per-band weight table rather than prefix sums, because a fractional edge
  cell's power sits at the midpoint of the *overlap* and weighting it by the bin centre can place a
  centroid outside its own band. One-sided bins are folded, so a full-scale sine inside a band reads its
  own mean square. Both axes are analysed: an anti-phase bass note would vanish from a Mid-only
  spectrum, and that is exactly the programme the first half is shouting about — so the peak band
  reports its own side fraction, i.e. whether the dominant note is cuttable.
  The background is the median **density** of the non-peak bands: semitone bands widen with frequency,
  so under a flat spectrum the top band of the range already holds 5.02 dB more energy than the median
  with no note present. The argmax is published twice (energy and density) and the dominance ratio is
  deliberately not a stored field — its denominator is exactly zero for a tone in digital silence.

Law 8a throughout (bit-identical under arbitrary re-slicing, stronger than law 11(a)): one integer
clock, samples outside and channels inside, `maxBlock` sizing nothing, absolute frame and block grids,
the denormal flush on `core::StateGrid` rather than at the end of `process()`, and no transformed tail —
`finish()` invents no frame and names `tailUncoveredSamples()` instead. Capacity exhaustion is data: the
block series keeps a prefix and says so while every integral, histogram and extremum keeps counting.
`fcore_measure lowend` prints the report as raw IEEE-754 bit patterns.

The suite nulls against oracles computed outside the object — the LR4 prewarped transfer function
predicts the settled side fraction at 120/180/240/480/1000 Hz, a direct O(N²) DFT with an independently
written band integration nulls the fold, and a full-scale sine pins the absolute calibration — because a
re-slicing test compares the implementation with itself and cannot see a deterministically wrong
schedule or a constant calibration error.

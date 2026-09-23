### Fixed — the ceiling density was the constant 1, and nothing asked it anything

`PeakExcursions::ceilingDensity()` and `ceilingDensityAbove()` shipped in v0.45.0 returning **exactly 1 for
every programme**. The histogram's bins run *downward* in level — `b = floor((12 dB − level) / 0.05 dB)` —
so the loudest populated bin has the *smallest* index and every quieter maximum sits at a larger one. The
loop walked `b = 0 … top`, which is the range *above* the loudest maximum and is empty by construction;
numerator and denominator were both the top bin. A consumer measured 1 on 102 files and said so. **There
was not one test for either accessor** — that is why a constant could ship as a measurement.

The loop now walks downward from the loudest bin, and both widths are clamped in `double` before becoming
an `int`: the wasm harness passes `1e9` dB for "every maximum", and `1e9 / 0.05` is not representable as an
`int`, so the conversion alone was undefined.

### Fixed — a flat top counted once per sample

The local-maximum test was `prev >= left && prev >= right`, under which **every** sample of a plateau is a
maximum. A clipped programme's L-sample flat top therefore contributed L counts at the ceiling, and the
density read the *duration* of the flatness as its crowding — on exactly the material the field exists to
describe. It is not a corner case: the reconstruction's group delay is 63.5 oversampled samples, half a
sample, so the crest of an isolated impulse falls *between* two bit-identical samples and even a lone spike
counted twice. The rise is now strict and the fall is not, so a plateau counts once, at its first sample.
`ceilingMaxima()` and both densities move for every programme.

### Fixed — a refused request answered with a number in range

A density lives in `[0, 1]`, so `0.0` is a legitimate reading — "nothing sits near the loudest maximum". A
bad argument, or a call before any measurement, returned `0.0` too and could not be told apart from one.
Both now return **−1.0**, as `LowEnd::sideFractionBelow` already does.

### Fixed — a crest figure that belonged to no fixture

"A two-sample click reads 10.8 kHz" in the K10 header and release note was carried over from elsewhere. A
click's excess over the ceiling grows with its amplitude while its width does not, and the estimator reads
`e/A`, so **there is no single click frequency**: it is **8.40 kHz at amplitude 0.85** and **11.40 kHz at
1.0**. Both are now printed by the suite that cites them, as is the 60 Hz sine's 59.3 Hz — the header said
59.4.

### Fixed — and the same half-convention one function over

`fc_probe_lowend_side_fraction_below` returned `0.0` when no measurement had been run, while the comment
directly above it declared the `-1.0` rule and the core already obeyed it for a refused frequency. A caller
that trusted the note read an un-run probe as a perfectly mono low end. It now returns `-1.0` there too.

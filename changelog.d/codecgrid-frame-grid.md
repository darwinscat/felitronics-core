### codecgrid — a new module: the frame grid a lossy codec leaves in decoded PCM

`felitronics::codecgrid` is an offline measurement. A decoder's output, analysed again with the codec's own transform
from the same sample, gives the quantiser's zeros back; one sample off, it does not. `CodecGridDetector` takes eight
stretches of a programme and reads each with seven hypotheses — the MPEG-1 Layer III hybrid filterbank, the AAC MDCT
with its sine and its Kaiser-Bessel-derived window, each at the programme's rate and at the other of 44.1 / 48 kHz
(the stretch is taken back with a zero-phase resampler), and CELT's MDCT at 48 kHz — at every frame offset. It reports
which grid stands out, how far, in how many stretches and on which phase: `Confirmed`, `InPlaces`, `SeveralGrids`,
`NoGrid` or `NotExamined`, with the family, the codec's rate, the hop, the phase counted from sample 0 of the
programme, every stretch's readings, the curve of the best reading and the zero share of its 32 bands with a map of
it per frame.

It reports a grid, not a codec: material that went through the same transform with coefficients zeroed by something
else reads the same. `NoGrid` means no grid of the transforms tried, on whole-sample offsets, in the stretches
examined — never "lossless". Long blocks only; HE-AAC, whose core runs at half the rate, has no hypothesis of its own.

A grid is `Confirmed` on either of two grounds. The rule (`GridRule`: score, depth in decibels, breadth over the 32
cells of four signals by eight band groups, uniqueness against the best offset more than two samples away) finds it in
two stretches on exactly the same phase; or the best offset of one hypothesis lands within a sample of one phase in at
least four of the eight stretches and the rule finds the grid on that phase in at least one of them. One stretch found
and nothing to confirm it is `InPlaces`.

A score counts the frames of a stretch as independent readings, and the frames of a programme that repeats on a whole
number of frames are not: sample-locked drums at 125 bpm put the same sixteenth note every six frames of 20 ms, the
noise of the curve is then that of six frames instead of a hundred, its largest offset scores like a grid, and it sits
on the same phase in every stretch. `GridScan::frameRepeat` compares the frames of a stretch with each other and
`repeatFactor` says what a score is worth among the frames that differ; a reading the rule accepts only at face value
is `repeated` and confirms nothing. The suite's rhythm — one noise burst every six frames, never coded — is found by
the rule at face value with scores of 17.6, 9.9 and 12.8 and is worth 4.5 at most. The rule's numbers, the stretch length of 2 s, the eight stretches at tenths of the programme and the phase tolerance are
parameters whose defaults were set by a measurement on real encoders that is not part of this repository; the suites
here prove the mechanism on programmes made to carry a grid — a transform, the smallest 30 % of every frame's
coefficients set to zero, the inverse — with no codec and no audio file.

The pieces are public: `Mp3Hybrid` (the polyphase analysis and the 36-point MDCT folded to an 18 x 18 DCT-IV, with the
alias butterflies), `Mdct` over a `ComplexFftBackend`, `MixedRadixFft` (the scalar reference for any size 2^a 3^b 5^c —
CELT's frame needs 480 points), the windows, `LevelAccumulator` and `readCurve`, `BackResampler`, `GridScan`. The
transforms run in single precision from tables designed in double through `core::det`. Against the standard's formulas
written out directly the MP3 analysis nulls 129 dB under the peak and the hybrid lines 134 dB; the MDCT nulls 134 dB
under the peak against its 2M x M sum for the three frames in use; `tests/crossnull` nulls the curves, the readings and
the resampler against numpy and scipy.

The form is the offline measurers': `storageFor()` publishes what `prepare()` takes and it is taken byte for byte
(10 029 524 bytes at 44.1 kHz and 10 591 124 at 48 kHz with the default stretches); `prepare()` is the only
allocation; `process()` keeps a copy of the stretches, takes any split of the programme and gives the same report for
every split; a request it cannot honour in full is refused whole; `finishStep()` does a bounded piece of the analysis
per call and `finish()` drives the same steps. Parameters it could accept and not honour — a phase tolerance wide enough to call any two offsets one phase, more
agreeing stretches than there are, a NaN anywhere — are refused. A rate other than 44.1 or 48 kHz, a programme shorter than one stretch
and a silent one are accepted and reported `NotExamined` with the reason. Non-finite samples are read as zero and
counted.

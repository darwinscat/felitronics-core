<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### analysis — `ProgrammeReport`: one offline pass over a delivered programme, and what it refuses to say

A new offline analyzer, `felitronics::analysis::ProgrammeReport`, shaped like `analysis::ClipDetector`
(`setParams` / `prepare` / `process` / `finish` / `reset`, `Storage` + `storageFor()` published before a byte
is allocated). One pass produces: DC offset per channel; leading and trailing silence, both at a named
threshold and as threshold-free exact digital silence; the tail — the energy of the last `tailWindowMs` of
programme relative to the programme's own mean square, plus the last sample and the last sample before the
trailing silence; the infra-low energy fraction through `eq::Crossover2`; the stereo relations (L/R balance,
bit- and value-identity counts, the integral correlation via `analysis::StereoSums`, side/mid energy); and
the dynamics — integrated loudness, reference true peak (`analysis::ReferenceTruePeakMeter`, 4× / 32 taps,
drained), PLR, EBU Tech 3342 LRA, the short-term loudness percentiles and the crest factor. Also
`fcore_measure report`, which prints every scalar as a raw IEEE-754 bit pattern with its validity and
reason, through the report's own field visitor.

**Every mean is taken over the PROGRAMME SPAN** — the first to the last frame whose loudest finite present
channel exceeds `silenceThresholdDb` — and the tail window ends at the last such frame rather than at the
last sample of the file. The property that buys: padding a master does not move its measurements. The DC,
sample peak, RMS, crest factor, infra-low fraction, programme mean square, tail ratio, stereo relations and
reference true peak of the same music padded with five seconds of silence at each end come back
bit-identical, pinned as an allow-list so a new padding-sensitive field cannot widen the claim quietly.
Measured over the whole file instead, every one of them moves for a reason that has nothing to do with the
music: a 60 s programme padded to 70 s reads 0.67 dB lower, and the tail window lands entirely inside the
padding and reads 0, so a truncated fade reports as a perfect one.

The span is thresholded rather than keyed on exact zero, which was the first design and is a cliff exactly
where masters live: a 24-bit dithered file — the dominant delivery format — has no exact zero anywhere, so
its five "silent" seconds are ±2e-7 of dither, an exact-zero span swallows them, the RMS is diluted by
5.4 dB and a programme cut mid-fade reads a tail ratio of 9.1e-12 instead of 1.0. What padding does still
move is named in the header: `lastSample`, which is by definition the last sample of the file, and the
loudness family, because BS.1770 anchors its gating blocks and 3 s windows at the start of the stream —
anchoring the sub-meters at the first active frame would make them invariant and put this report's loudness
at odds with `fcore_measure lufs` and every other EBU tool on the same file, which is the worse trade.

**"Cannot say" is a result.** A mono file does not get a correlation of +1.0; a programme shorter than the
tail window does not get a ratio of exactly 1.0; a programme with fewer than two gated short-term
observations does not get an LRA of 0.0, which is also the honest answer for a constant tone. Each is
`valid = false` with a reason and a canonical `+0.0`, never a NaN. Two shipped primitives are wrapped
because their own defaults are verdicts rather than measurements: `StereoSums::correlation()` answers +1.0
below its 1e-12 denominator ("silence is neutral" — right for a meter, wrong for a report), and
`LoudnessMeter::integratedLufs()` answers −120.0 when nothing passed the gates. LRA is computed from this
class's own short-term series rather than from `loudnessRangeLu()`, whose 0.0 cannot be told from outside
the meter to mean "no range" or "fewer than two observations survived the gates".

**The infra-low fraction is a filter, not a band, and the header says so in the arithmetic.**
`eq::Crossover2` is LR4 — two cascaded Butterworth sections — so |LP|² = 1/(1+(f/fc)⁴)², which is −6 dB at
the crossover, **not** the 4th-order Butterworth 1/(1+(f/fc)⁸) the request was written against. With
fc = 30 Hz a pure 30 Hz tone reads 25 % rather than 50 %, a 40 Hz tone 5.78 % rather than 9 %, and flat
noise 0.104 % at 48 kHz. A reader told the wrong formula would see a few per cent "below 30 Hz" on a
bass-heavy master and conclude there was infrasound, where it is the fundamental at 40 Hz through the
filter's skirt.

**Law 8a, gated on the intermediates and not only on the result.** Every field, and every event of a
test-only frame trace, is bit-identical under any re-split of the stream into `process()` calls — checked
over 19 slicings, 8 `maxBlock` values and anchored random slicings, field by field through `std::bit_cast`.
The trace exists because the finished report is not a sufficient gate, and that is measured rather than
argued: of a 23-mutant stand, 22 die and **not one of them produces a single report-comparison failure**.
A mutant that closes the 10 ms sub-hop at the end of `process()` leaves the report bit-identical at all 19
slicings and all 8 `maxBlock` values; so does one that moves `eq::Crossover2::flushDenormals()` off
`core::StateGrid`. Beside the invariance sit closed-form oracles computed outside this repository (a sine's
crest factor is 20·log10√2 = 3.0103 dB; the LR4 magnitudes above; EBU Tech 3341's −23 dBFS tone) and an
independent whole-file reference program that takes its percentiles by sorting rather than from a histogram
and writes out the 3 s / 1 s cadence as literals rather than reading the header's constants — an oracle
that imports the subject's parameters moves with a mutation of them instead of opposing it. Two real
defects in this work were caught by that reference alone and by no invariance check.

Non-finite input is a hole: a canonical zero into every filter, counted per channel, and then excluded where
exclusion is exact (the sums stay valid) and fatal where state carries it forward (the infra-low fraction
and the whole loudness family refuse, with the reason). Nor is a maximum exclusion-exact — the sample
dropped for being non-finite may have been the largest — so with nothing finite the sample peak is refused
rather than published as a comfortable zero. Finite input is not enough either: `eq::Svf` updates its
integrators as `(float)(2·v − ic)`, and that intermediate can leave the float range even where the state is
not growing, so a stream of finite `FLT_MAX` drives the 30 Hz crossover non-finite at sample 849 exactly,
grid flush running (a 1.7e38 DC input does not overflow, so "the state doubles every sample" is the wrong
mechanism to quote). The crossover's overflow and the K-weighting's are counted separately, so a refusal
never names a filter the field does not go through. Nothing non-finite reaches an accumulator, a published
field or the trace. Capacity overflow never refuses a call or stops a counter; it invalidates exactly the
fields it damages.

475 checks, green under Apple clang 21 / libc++ / arm64 and gcc 14.2 / libstdc++ / x86-64, including the
strict header-hygiene gate (`-Wconversion -Wfloat-equal -Werror`) on both. `fcore_measure report` on a
stereo fixture is bit-identical across those two rows — every count and all 31 scalars — which is measured
here rather than promised: the header does not claim cross-platform bit-identity, since `log10`, `pow` and
`sqrt` are in the derivation.

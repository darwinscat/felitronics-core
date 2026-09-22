### analysis · probe — BandCrest: the crest a band lost between a source and its master

The micro-dynamics damage counter an automatic mastering mode needs. Macro damage is LRA loss and already
exists; timbre is a band energy shift and does not yet; this is the third. Run it on the source, run it on the
delivered master, and the damage is the crest a band lost, block by block.

**What a band crest can and cannot mean** is in the header, because the honest limits are part of the
instrument: it is the peak-to-RMS of a signal seen through a NAMED filter, never a statement about "the content
between 2 and 6 kHz". The band is a dome and an LR4 skirt is only 24 dB/octave down, so a band's peak can be a
neighbour's transient leaking in. In the DIFFERENCE of two runs through the same filter that leakage largely
cancels; in one run's absolute number it does not.

#### Four decisions, each where the obvious choice is wrong

**Reconstructed peak, not sample peak — and the reason is not the sample rate.** Even with both runs on the
same grid a sample peak is biased in ONE direction: the input carries inter-sample peaks and the output has
been through a true-peak limiter whose whole job is to remove them. A sample-peak crest would under-read the
input's crest and therefore UNDER-REPORT the loss, which is the one number this exists to produce. The
reconstruction is the certifying instrument's own 4x/32, so the full-band peak can be nulled against the
certificate rather than believed.

**One reconstruction, then the split at 4x** — cheaper than five reconstructions, and it makes the band
magnitudes nearly rate-independent as a side effect, since an LR4 corner is prewarped by `tan(pi*f/fs)` and at
the base rate the skirt moves with the sample rate.

**The bands are independent filters of the input, not a reconstructing tree.** A tree is cheaper and sums back
to the input, but pays for that with allpass compensation on every band — and an allpass changes the waveform,
which is to say the peak, which is to say the very quantity being measured. Nothing here needs the bands to
sum, so nothing here pays for summing. They do not reconstruct the input and are not meant to.

**The activity mask is built from the SOURCE only.** This is the decision that can quietly invalidate
everything else: a mask taken from the output, or from each side separately, can drop exactly the blocks the
chain damaged most — a loud block the chain crushed can fall below a gate its input passed, and the damage
leaves the population with it. The gate is a conjunction, because each half fails in its own direction: an
absolute programme floor (a relative loudness gate would throw away the quiet passages where a compressor
works) AND a band share floor (without it a band with no content contributes the crest of its own noise floor,
and a source that is digital black in a band against a master carrying dither reads as several dB of "damage").

#### Four numbers, not one — the obvious one is blind to the instrument doing the damage

A limiter works on the loudest few percent of blocks. The suite's own fixture has **p95 reading 0.000 dB while
CVaR95 reads 2.214 and 8 of 297 blocks lost more than 3 dB**: a cost function calibrated on p95 would call that
master undamaged. `cvar95Db` — the mean of the worst 5 % — is the number to read; `p95Db` stays beside it
because its disagreement with CVaR95 is itself the signal that damage is concentrated; the exceedance counts
say the same as a count; and `p5Db` is the negative side, where a limiter's pumping shows (it lowers a band's
RMS around a kick while the peak elsewhere is untouched, so the crest RISES and the positive part hides it).

`peakShiftDb` and `levelShiftDb` ride along so a crest loss with no peak loss — a floor RISE, a clipper putting
a steady harmonic into a band that held only a click — can be told apart from a flattening.

The requester asked for the check in the other direction too: a fixture where p95 sees damage and CVaR95 does
not must not exist. It cannot, by construction — CVaR95 averages the worst 5 %, every member of which is at or
above the 95th percentile — but "cannot by construction" is exactly the kind of statement that is false at the
edges of a rounding rule, so it is checked over 400 shapes including populations of one, all-zero, all-equal
and a single outlier.

#### Every oracle is outside the object, and two of them caught defects before any test was written

The full-band peak is **bit-identical** to `ReferenceTruePeakMeter` at 44.1, 48 and 96 kHz — the same topology,
so a tolerance would hide a difference rather than allow one. The block grid equals `LoudnessMeter`'s gating
block count at all three. The interpolator's delay is measured from an impulse: symmetric about oversampled
index 319.5 with a worst difference of exactly 0, which is where the 63.5 in the header comes from — a HALF
sample, named rather than rounded away.

Two independent hop clocks, counted by different code over the same time, **disagreed 31 against 30**: the
interpolator's drain was being taken whole, so a hop of silence past the end of the programme closed as
programme. And the RT check passed on clang and failed 25 times on gcc — `alloc::count.load()` was one argument
to `ok` and an allocating `std::string` was another, and argument evaluation order is unspecified.

#### The surface

`fc_probe_crest_run[_with](slot, ...)`, `_scalars`, `_blocks`, `_loss`, `_storage_bytes[_with]`, plus
`fcore_measure crest [--against other.f32]`. **Two slots, not one**, because the loss is arithmetic that must
not be reimplemented across the ABI: formed from linear cells it is one logarithm and invariant to the gain a
master has over its source, and reassembled in JavaScript from two dB columns it would be neither. The
per-band accepted-block counts are published rather than left to be counted from the block table — a count a
consumer derives itself is a second definition of "active" waiting to drift.

An unrecognised `--` token is a refusal on both roads from the start. Native and wasm are byte-identical on the
default road, on a parameterised one and on a run carrying the loss rows; their refusal sets agree on twelve
argument lists, ten of them refusals.

#### What two review rounds found, and what that says about the tests

A code review found seven things and a falsification round — whose mandate was to make this analyzer report a
wrong number rather than to read it — found nine more. Both are worth listing, because six of the sixteen were
cases where a TEST was narrower than the claim it stood under.

**The full band's crest was a hybrid and wrong near Nyquist.** Its peak is the reconstruction; its mean square
came from the oversampled stream, which the 32-tap interpolator droops above about 0.375·fs. A pure sine has a
crest of exactly 3.0103 dB — an oracle that owes this code nothing — and it read **9.03 dB at 0.45·fs and
23.8 dB at 0.49·fs**. The per-band crests were never affected, because their peak and their mean square come
from the same stream; only that one number was wrong, and it fed the absolute gate and the published programme
level with it. The full band's mean square is the base-rate one now and the share gate keeps the oversampled
denominator, because a ratio across two domains is a share of nothing.

**The refusal set differed by platform.** Six sibling analyzers carry a 768 kHz ceiling; this one had none, so
`llround` decided — and it saturates on arm64 while it wraps on x86-64 glibc. Measured on both rows:
`fs = 1e308` was refused on one and accepted on the other with a 10-sample hop, answering a measurement that
looked entirely valid. The bounds are comparisons made before any conversion now, and the two rows agree.

**The filters were `SystemMath`**, where every sibling is deterministic and the math-policy suite asserts it of
each — `std::tan` disagrees between libms and moved 12 block rows of 3752 on a 6.3-minute mix. **The hop was
`llround(fs·hopMs/1000)`** where the loudness meter builds it from sub-hops, so 22050 Hz gave 2205 samples
against 2210 — and the test that claimed the identity ran at 44.1, 48 and 96 kHz, the three rates where the two
agree. **A refused `prepare()` answered the previous run.** **The comparator compared anything**: the same audio
with the master on a 50 ms hop gave CVaR95 7.8 dB with `valid` true. **No evidence was reported as a lag**: a
muted master read "misaligned by 8 blocks". **The non-finite count meant interpolator outputs**, 128 for one
NaN, and its location depended on the call sizes. **A master scaled by 2⁻⁵⁰** — float subnormals, strictly
positive — was compared as real. **A channel that stopped was not cleared**, so a returning one spliced a
32-sample-old ring onto new audio: 40 dB of crest on a block. **And the budget was a tautology**: `bytes() > 0`
passed while the demand omitted the interpolators themselves, 168 B a channel at every length.

**The bed placed every transient on the hop grid.** At t = k·0.5 s, which is a hop boundary at every rate this
suite uses — the feature sat exactly on the grid the measurement is made of, in the fixture the alignment tests
are built on. That is this repository's own rule 2 inverted. Moved off it by a third of a hop, which turned the
CVaR95 control red, because that control's threshold was a number fitted to the fixture of the day rather than
a statement of the mechanism. It states the mechanism now.

Also measured and then gated: gain invariance is **exact only where the scaling is exact** — ×0.5 is, ×0.37 is
not, and the 1.1e-5 dB that appears there is the fixture rounding rather than the comparison drifting, so the
claim is now the two claims it always was. And the filter bank runs the twelve evaluations it reads instead of
the twenty it computed: **2.84–3.12 s → 2.04–2.05 s** on a 4-minute stereo programme, bit-identical — but only
after a diff said it was not, because the first version passed `kQ` cast to float where `Crossover2` passes a
double.

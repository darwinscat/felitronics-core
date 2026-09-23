### Added — K13, a peak clipper, and it lives INSIDE the limiter's oversampling island

A true peak clipper that shaves the tallest, shortest excursions before the limiter's detector ever sees
them. `TruePeakLimiterParams` gains `peakClip`, `overCeilingDb` and `kneeDb`; the ABI gains `peakClipper`,
`peakClipperOverCeilingDb` and `peakClipperKneeDb` (v11).

**It is not a chain stage, and that decision is the whole design.** The limiter already upsamples, acts on
the oversampled grid and downsamples — the structure its header calls "Option B, the only one that can
bound inter-sample peaks at all" — and that island carries the *audio*, not just the gain decision. A clip
placed inside it costs **nothing**: no second round trip, no added latency, no `DryAligner`, and not one
field of the PDC moves. A separate stage before the limiter would have cost a third oversampler round
trip, 63 samples, and −15.55 dB at 20 kHz (44.1 kHz, the oversampler's own figure) — paid on every
programme, including the ones where the clipper never engages.

**Two different ways of not clipping, and both are bit-exact.** `peakClipper = 0` never enters the code;
an `overCeilingDb` so high that nothing reaches it enters and finds nothing to do. They are not the same
setting — the first cannot start clipping when the ceiling moves, the second can — and a page that means
"off" should say so with the flag. The whole existing suite is green with the clipper off, which is the
anchor: an old parameter set renders exactly as it did.

**The threshold rides the ceiling.** `overCeilingDb` is an *offset above* `ceilingDbTp`, in [0, 12], not an
absolute level — so sweeping the ceiling cannot leave the clipper under the thing it guards. The absolute
level it lands on, after both clamps, is published as `fc_master_resolved::peakClipperThresholdDbTp`.

**The knee is C1 at both joins**, quadratic in the linear domain, monotone, gain never above 1, and an
exact hard clip at `kneeDb = 0`. A hard clip leaves a first-derivative corner whose harmonics fall off like
1/n²; this leaves a second-derivative corner and 1/n³. That is real alias rejection while the excess is
inside the knee and nothing once the plateau dominates, which is why the range stops at 1 dB.

**Channel-linked**, one gain for all channels, like everything else here: a per-channel clamp satisfies the
same bound and moves the stereo image on exactly the samples a listener notices.

**What it promises and what it does not.** Promised, and it is algebra: every sample on the limiter's own
4× grid leaves the clip at or under the published level. **Not** promised: that the delivered file's true
peak is under it. A hard-clipped sine comes back from the downsampler as its own fundamental *above* the
clip level — the limit is **4/π = +2.10 dB**, the square wave's first Fourier coefficient, and no
oversampling factor touches it. The same gap the limiter's own header already describes; the clipper
inherits it rather than adding a new one. `truePeakDbTp` is what shipped.

**`limiterMaxReconstructedPeakDb` is unchanged and still means the peak that ARRIVED**, before the clip —
the tap it is read from is written one line above the clip, and a test pins it. How much the clipper took
off is the difference between that and the output.

The measurement gains six fields on the oversampled grid, in K10's terms: `peakClipReductionMaxDb`,
`peakClipReductionP95Db` (a quantile over *clipped samples*, 0.1 dB resolution, not a window quantile),
`peakClipOccupancy` (−1.0, never 0.0, when nothing was judged), `peakClipRuns`,
`peakClipRunSamplesTotal`, `peakClipLongestRunSamples`.

Cost: the limiter object grows 1 KiB for the reduction histogram the quantile is read from. 0.1 dB rather
than 0.05 for that reason — at 0.05 it would have been 2 KiB against an object of 744 bytes.

### Fixed — the ordering the design was written against, pinned before it could ship

The clip must sit **after** the two lines that record what arrived and **before** the detector. One line too
early and `limiterMaxReconstructedPeakDb` reports the *clipped* value: on a programme arriving at
−0.10 dBTP with the clip at −6, the field would have read exactly −6 — always the clip level, a constant
dressed as a measurement, plausible and under the ceiling and invisible to any range check. The planted
control is in the suite.

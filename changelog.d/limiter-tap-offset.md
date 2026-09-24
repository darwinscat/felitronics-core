### Fixed — `limiterTapOffset` ran 0.875 frames early, and the test could not see it

`resolved().limiterTapOffset` published **31** where the truth is **31.875** frames. A consumer cropping a
statistic by it cropped a frame short; the crop now moves by one frame, which is the one visible change.

**The up leg is not half the round trip.** The limiter's trace is written where the gain is *decided* — on
the oversampled copy — so the tap lags by the up leg alone, and for the Kaiser polyphase that is
`(N − 1) / (2F)` with `N = factor × tapsPerPhase`: 255/8 = **31.875** frames at 4×/64. The expression took
half the *reported* round trip instead (63/2 = 31.5) and then truncated it to 31, so the integer division
was only part of the error.

**Measured, not argued**, because two earlier readings disagreed. With the limiter neutralised and an
impulse into the chain, the tap's response is exactly symmetric about oversampled index *x*.5 — the samples
either side are bit-identical and so are the pairs around them — and its energy centroid lands on
**31.8750** frames. The maximum *sample* sits at 31.75, a quarter frame early, because a half-sample group
delay puts the crest between two equal samples. Reading the argmax is what made this look like a simple
truncation, and it is why the test now reads the centroid and asserts that the two disagree: if they ever
agree, the half-sample delay has moved and the tolerance must be re-derived.

**The audio was never affected.** The limiter's own delay measures exactly 111 samples and its output
impulse response is exactly symmetric, so the chain's PDC was right all along — only the published tap
coordinate was wrong.

The test's tolerance was **1.0 frame** and therefore could not tell 32 from the 31 that shipped. It is now
**0.5**, which is what a contract in whole frames can mean at all; the residual 0.125 is the distance from
the truth to the nearest frame and nothing about the measurement.

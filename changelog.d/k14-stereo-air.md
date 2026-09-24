### Added — K14 "stereo air": a high shelf on Side, inside mono-bass's own M/S island

`stereo::StereoAirParams` — `enabled`, `frequencyHz` and `gainDb` — carried by `stereo::MonoBass`, which
already opens the island. The ABI gains `fc_master_config::stereoAir`, `fc_master_params::{stereoAir,
stereoAirHz, stereoAirDb}` and `fc_master_resolved::{stereoAirHz, stereoAirDb}` (v12).

**It rides the existing island rather than opening a second one.** The M/S round trip is not the identity
in float — `(m+s) + (m−s)` rounds each term separately — so a second one would perturb the programme twice
for nothing. The island is now entered when **either** tool is configured, both early exits ask about both,
and each tool keeps its **own** latch: resetting the crossover because the air toggled would click the
bass, and resetting the shelf because the bass settled would click the top. A mono chain is refused with
`stereoAir` exactly as it is with `monoBass` — a stereo tool silently doing nothing is the class this chain
already closed.

**`gainDb` is the PLATEAU.** The shelf reaches **half** of it at `frequencyHz` and the rest above: measured
on the tree's own designer at 48 kHz with a 6 kHz corner and a +3 dB request, 1 kHz sees +0.0024 dB, the
corner **+1.4983**, 10 kHz +2.6479 and 20 kHz +2.9801. The band *energy* moves by **×1.9135** for that
request, not by the plateau's ×1.9953, because the measurement band is an LR4 high-pass — a weighting, not
a wall. `frequencyHz` is clamped to `[3000, min(12000, 0.45·fs)]` and read back through `resolved`.

**Zero is a branch, not a unity filter.** `matched::highShelfDb` *substitutes* `g = 1.00001` for a zero
request — a 0 dB shelf is designed as a +8.7e-5 dB one — and even forced identity coefficients are not
transparent: the biquad computes `1.0*x + 0.0` in double and turns **−0.0f into +0.0f**, measured, which is
the defect that disqualified the clipper's `mix = 0` bypass. So the stage is skipped by a predicate on the
**smoothed and clamped** value, in the shape `StereoWidth` already uses, including its check that the ramp
has settled. The flag and a 0 dB request take the same path, and both are bit-identical to the tool being
absent — a planted control that runs the shelf at zero differs on **47767 of 48000** samples.

### What it does to the sound, and what it does not

- **The mono fold does not change.** `(l+r)/2` is Mid, and Mid is untouched — for any Side processing,
  phase included. What grows is the **gap**: on uncorrelated highs a +3 dB plateau takes the stereo-to-mono
  drop from −3.01 dB to −4.76 while the mono level stays exactly where it was.
- **On anti-phase highs every width number is blind.** Width reads **1.000 before and 1.000 after** while
  the Side energy grows by the full band integral, because the fold was already empty. That is why the
  measurement publishes **three energies** (`airMidEnergy`, `airSideEnergyBefore`, `airSideEnergyAfter`)
  and not only the fraction; the two `airWidth*` fields use the page's own amplitude convention
  `√S/(√M+√S)` so they sit on the same scale as its broadband meter, and are **−1.0** — never 0.0 — when
  there was nothing to judge.
- **It bleeds a hard-panned top into the other channel.** With content only in L, the output R comes up to
  about **−15.4 dB** of L with **inverted** polarity above the corner, because `r' = m − H·s`. No L/R shelf
  does this; it is a property of acting on Side, not a defect.
- The mono **sum** is preserved to within **one ulp**, not bit-exactly — the same rounding the island's own
  header already documents. Where Side is exactly zero (`l == r`) the reconstruction *is* bit-exact.

### Fixed — a topology flag that reached neither the core nor the hand mirror

`config.stereoAir` was copied into `MasteringChainConfig` by nothing, and into `fcore_master`'s independent
C++ mirror by nothing either: a page could have set it and got a silent no-op — the exact failure the flag
exists to prevent. Two different gates caught the two halves, the domain table and the ABI selftest.

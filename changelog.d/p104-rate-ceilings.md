<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — the search and the shapes get the ceilings their neighbours already had

**What each entry now promises about the rate, from above:**

- **`mastering::TargetLoudnessSolver`** measures up to 3 MHz — `kMaxSampleRate`, which is
  `MasteringChain::kMaxSampleRate`: the search measures what a chain renders, at the chain's rate. `prepare()`
  refuses a higher rate and disarms, and `solveBytes()` / `measureRangeBytes()` answer 0 there. It took any
  finite rate over the floor: at 1e300 Hz `prepare()` said yes while `solveBytes()` said 0, and one ulp over
  3 MHz it prepared for a chain that cannot exist.
- **`fcore::ShapeProbe` and `fc_probe_shapes_run`** take 8000 Hz to 768 kHz — `fcore::Probe`'s range. They
  kept no ceiling and drew a 1 MHz file that every other run entry of the probe ABI refuses. `fcore_measure
  waveform` / `stereo` refuse above 768 kHz too, and name the range.

Nothing changes at or under either ceiling: a verdict oracle over 3817 rates and 39 entries finds differences
only in the search's and the shapes' entries, only above their ceilings, and only from an acceptance to a
refusal.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### core · analysis · mastering · tools — one sample-rate floor, 8000 Hz, and a rate in kilohertz is refused

**What each entry now promises about the rate** — one number, `felitronics::core::kMinSampleRate = 8000.0`,
compared as `sampleRate >= kMinSampleRate` (8000 is a rate; NaN is not):

- **`mastering::TargetLoudnessSolver`** measures at 8000 Hz and up. `prepare()` refuses a lower rate and
  disarms; `solveBytes()` and `measureRangeBytes()` answer 0 exactly there. It took any finite rate > 0.
- **`mastering::MasteringChain`** is built at 8000 Hz up to its own 3 MHz, on every topology — the floor is the
  chain's, not a stage's. Before, what refused a low rate was whichever stage happened to be on (the limiter
  above 50 Hz, the EQ above 20.4 Hz), so a chain without both took 1e-305 Hz. `admits`, `prepareBytes`,
  `reprepareBytes` and `mastering::createBytes` say the same.
- **`core::DeliveryResampler::plan`** plans integer rates from 8000 Hz (its literal floor was 1000), so
  `mastering::DeliveryConverter` and `DeliveredMastering` refuse a SOURCE or a DELIVERY rate below 8000 — in
  `prepare()`, in every budget, and in `deliveredFrames()`, which answers -1 there. This is the line that covers a delivering handle's source rate: its chain
  runs at the delivery rate and never sees the source.
- **The mastering C ABI:** `fc_master_create` and `fc_master_need_create` answer `FC_ERR_REFUSED_BY_CORE` for
  `sampleRate` (or a non-zero `deliveryRate`) below 8000, with nothing allocated. No ABI version change: the
  surface did not grow, and the refusal comes from the core the facade already asks.
- **The probe ABI and `fcore_measure`:** `fcore::Probe`, `fcore::ShapeProbe` and the six analyzers behind
  `fc_probe_*` — `ClipDetector`, `ProgrammeReport`, `HumDetector`, `LowEnd`, `SourceForensics`,
  `BandBursts` — take 8000 Hz and up (to 768 kHz; the shapes have no ceiling). The six analyzers and the probe had a floor of 1000 Hz, each
  in its own copy, and the shapes had none; each class keeps its `kMinSampleRate` name, and every one is
  now the core's constant. `fc_probe_<mode>_storage_bytes` answers 0 below the floor. The ceilings did not
  move, and the shapes still keep none (a rate Probe refuses above 768 kHz is still drawn). `fcore_measure`'s
  `correlation` and `needle` modes do not read the rate and still take any finite positive one; every other mode
  refuses below 8000 and names the range.
  A refused `fcore::Probe` or `ShapeProbe` now reads like a fresh one — its meters and its peak and stereo parts
  are replaced by new ones, and their memory is released — instead of serving the previous file (the rates that
  used to re-prepare them are refusals now), and
  `DeliveredMastering::sourceRate()` / `deliveryRate()` read 0 after a refused prepare, as the chain's and the
  search's rates already did.

**What starts to be refused, and why that is a finding and not a loss.** Every call below 8000 Hz that was
accepted before — the probe and the analyzers between 1000 and 7999 Hz; the chain, the search and
`fc_master_create` anywhere below 8000 their stages let through (60, 88.2, 96, 192 Hz on the default
topology; 44.1 with the limiter off); delivering handles from or to 1000…7999 Hz; the shapes at any positive
rate their decimation could hold. Measured on the tree before this change:

- the BS.1770 K-weighting shelf is designed at 1681.97 Hz, so below 3364 Hz it is past Nyquist — aliased
  everywhere there, and **unstable** wherever the bilinear tan() comes out negative (1682–3364 Hz, 841–1121 Hz,
  …): `fcore_measure lufs 3300 2 fixture.f32` (the CI fixture) printed **+3048.86 LUFS**; a 0 dBFS 400 Hz sine
  reads −3.72 LUFS at 48 kHz and +3043 LUFS at 3300 Hz;
- a search at 3363 Hz answered TargetUnreachable **at its −60 dB gain rail**, chasing +2448 LUFS;
- a search at **88.2 Hz** — 88.2 kHz spelled in kilohertz — answered **Solved at −14.0 LUFS**, exit 0, with a
  whole mastering chain running a thousand times too slow; 384 and 768 did the same; `fcore_master render 96`
  wrote its file with no status at all.

8000 Hz is the lowest standard audio rate, so no real programme sits below it, and it is clear of the shelf's
edge (its pole radius is 0.43 at 8 kHz, 0.99997 at 3364 Hz). A `static_assert` in `KWeightingFilter.h` keeps
the floor above twice the shelf. **Not changed:** `analysis::LoudnessMeter` and `KWeightingFilter` still take any
rate, and the meter still reads a rate ≤ 0 as 48 kHz; every caller that takes a rate from outside is floored
above them, and the meter's own contract is a separate decision. A page that builds a Web Audio context below
8 kHz on purpose is refused too.

**Tests that ran below the floor moved above it, with their fixtures kept:** `HumDetectorTests` from 3 kHz to
12 kHz (order 15 — the same 0.366 Hz bin; the noise is the 3 kHz draws, interpolated, and every measured
prominence moved by at most 2e-4 dB), `LowEndTests` from 6 kHz to 12 kHz (order +1, every sample count ×2, the
noise of its band-statistics rows interpolated, the kick's click scaled to its per-bin power; the fixtures that
are only compared with an oracle computed from the same samples draw fresh noise), `ClipDetectorTests`'
short-stream null and pending-queue rows and `ClipsExposureTests` from 1 kHz (W = 20) to 8 kHz (W = 160), their
quiet gaps kept at the same number of windows. P41 F1's solve at 1e-305 Hz is unreachable
now; its property — a meter sized in samples — is pinned on `LoudnessMeter::prepareForSamples`, which still
takes such a rate. **CI:** the clips gate refuses 7999, one ulp under 8000, 44.1 and 1000, and measures 8000;
the price table's rates are 999, 1000, 7999, 8000, 11025, 12000, 16000 … and its floors are 475 rows / 284
priced (425 / 280 before — the 1000 Hz rows are refusals now, and the 2000 Hz rows left the grid).

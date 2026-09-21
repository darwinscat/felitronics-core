<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — the EQ points of the mastering chain are DYNAMIC: `eqBands[].dyn` renders

**No ABI change.** `FC_MASTER_ABI_VERSION` stays 8 and no struct moves. The six `fc_eq_dyn` fields were
already carried, clamped and handed to the core; what they were not was rendered. `eq::EqBand` applies a
gain delta a PRODUCER pushes into its per-lane seam, and the chain's EQ stage had none — so a point with
`dyn.on = 1` was, measurably, the same samples as one with it off.

**`mastering::MasteringChain` now drives one `dynamiceq::LaneDynamics` per EQ band.** The detector key is the
EQ section's INPUT, captured through `eq::EqEngine::captureSectionInput` before any band moves the signal, so
one point's moving delta cannot modulate a later point's detector at an overlapping frequency. Each point's
five placement lanes get their own probe, programme estimate and delta from one shared setting, which is the
composition `dynamiceq` was written for; nothing new is computed here.

**The default render did not move.** With every point unarmed — `dyn.on = 0`, which is the default — the EQ
stage is `eq::EqEngine::process` over the same quantum, sample for sample, and a producer whose point is not
armed runs its band over the WHOLE quantum in one call rather than in control-rate chunks. Checked against
v0.39.0 through the C ABI on twenty seconds of real programme with three points, mono-bass, compressor,
clipper, limiter and 24-bit dither: byte for byte identical, including a run with every `dyn` field set to a
far end while `dyn.on` stayed 0. `dyn.rangeDb = 0` is no dynamics whatever `dyn.on` says, and renders
identically too.

**It costs no latency.** The producer applies the delta derived from the PREVIOUS control chunk, so nothing is
read before it is written: `latencySamples()`, `compressorTapOffset` and `limiterTapOffset` are unchanged, and
an impulse leaves at the same index armed or not. Block invariance survives with a point armed — the
producer's 16-sample control grid restarts at every call, and the chain's call is always the internal quantum,
so the restart lands on the same absolute sample at every caller block size (measured bit-identical at 1, 7,
63, 256, 1021, 4096 and 65536).

**And the loudness search is undisturbed**, which is the property the EQ's place in the chain buys: the stage
sits ahead of the node the search moves, so the PRE-LIMITER tap is bit-identical at every drive with a point
armed, the limiter's worst gain reduction stays monotone in that drive, and `y(g, c) = 10^(c/20)·y(d, 0)`
holds to 4.9e-07 at a peak of 0.891 on a 3x3 grid — the same float-rounding order `LoudnessSolver.h` states
it at for a chain with no EQ at all.

**Memory is declared before it is asked for**, as every other stage's is: `Storage::dynBands` counts the
producers, `bytes()` carries them, and `prepareBytes()` is still what a fresh `prepare()` allocates byte for
byte on every row of the rate x width x topology matrix. One array of 24, held on the heap beside the engine
for the same reason — about 60 KiB on a chain that has an EQ, nothing on one that does not. `process()`
allocates nothing with points armed, arming and disarming included.

**The four `dyn` domain rows are pinned by behaviour now.** `FC_DOMAINS` carried `rangeDb`, `thrDb`, `atk` and
`rel` as clamps with no read-back — they were the last four rows of the table that nothing could pin, and the
reason was that the fields were inert. They now read back through a render, and `MasterDomainsTests` asserts
the opposite of what it used to: an armed point is NOT the same samples as `dyn.on = 0`. Pinning them needed
armings that make each bound visible — a threshold every sample is over for the `rangeDb` cap, an input gain
that moves the programme to meet each end of `thrDb`, and a threshold inside the programme's own band level
for the ballistics, since a reduction pinned at its cap never releases and a release knob then changes
nothing. Zero of 89 rows are now unpinned.

**`fcore_master` can drive it**: `band<i>.dyn.on`, `.dyn.range`, `.dyn.thr`, `.dyn.thrAuto`, `.dyn.atk`,
`.dyn.rel`.

**Not yet accepted by ear, and shipping that way on purpose.** Oleh listened to the A/B on 2026-09-21 — it
works, but a dynamic point is hard to hear on programme material, so this feature is lightly tested where
only ears can test it. That listening continues in TabbyEQ, where the same dynamics are still experimental;
until it concludes, treat the feature as one the measurements vouch for and the ear does not yet.

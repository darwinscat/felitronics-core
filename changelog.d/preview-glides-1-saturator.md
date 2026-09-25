### saturation — the Saturator's parameters glide, and the first write of a stream snaps

A write to `driveDb`, `bias`, `autoComp`, `mix` or `outputDb` used to apply at once, so a knob turned on a live
stream stepped the output: through the mastering chain on a -12 dBFS 227 Hz sine (K = 128), max|Δ²y| where the change
reaches the output was -18.3 dBFS for drive 3 -> 9 dB and -44.6 for a 3 -> 3.5 dB nudge, against -67.2 for the steady
saturated tone; the bias of the Asym curve 0 -> 0.3 read -11.1, trim 0 -> -3 dB -50.8, autoComp 0.5 -> 1 -48.6, mix 1
-> 0.5 -54.6. A continuous write now starts a linear ramp in PARAMETER space (`Saturator::kGlideMs` = 30 ms), stepped
on the 64-sample `core::StateGrid`: at each grid boundary the curve is designed at the ramp's next point with the
settled stage's own arithmetic, and inside the period the curve's coefficients (k, bias, tanh(k·bias) and the RAW
PEAK the normaliser is the reciprocal of; oversampled rate) and drive-compensation, mix and trim (base rate) are
interpolated per sample between the two designs — as `start + d·index`, product and sum in separate statements (law
10), so a coefficient is a function of the sample's place in its period and never of a call. The curve is divided by
the interpolated raw peak rather than multiplied by an interpolated normaliser: the normaliser goes like 1/k at a small
drive, and the review round measured a 0 -> 3 dB glide on a constant 0.2 peaking at 7.665 when it was interpolated
directly; a glide now stays inside the range the settled stage covers along the same path (pinned for every shape). Measured the same way: -63.8, -67.2, -63.6, -67.3, -65.2,
-65.3 — at or within 3.4 dB of the steady tone's own reading. 30 ms was measured against 10, 20 and 50: every length
from 20 ms is clean, 10 ms leaves the bias move at -59.7, and 30 ms matches the mastering chain's own ramps.

The first write of a stream SNAPS: every `setParams()` between `prepare()`/`reset()` and the first accepted call with
samples applies at once, exactly as before, so a stage whose parameters were set before its first sample renders the
bits it always rendered — the whole offline contract of the mastering chain above it (six chain topologies compared
by hash, identical). `reset()` lands a glide in flight. `shape` is a topology switch and snaps every parameter with
it; `dcBlockHz` lands at once — both as they always did. `setGlideMs (0)` turns the glide off and IS the pre-glide
stage, bit for bit: the frozen pre-change engine's mid-stream sweep (the hardening NULL's scenario 5) now pins exactly
that. A write that changes nothing continuous costs a comparison (it used to cost a full re-design, and `dcBlockHz` is
re-designed only when it moved), so a caller that re-sends its parameters every block pays nothing. The ramp accumulates in double, landing on the target itself. New
`WaveShaper::coeffs()` and `WaveShaper::shapeAt<S>()` expose the curve at explicit coefficients — additive, and
`processSample()` is untouched (the guitar core's TubeStage uses only that). `glideTicks()` and `isGliding()` read the
state back. The settled slice keeps its code verbatim and the glide runs in a separate function, so a stage that is
not gliding pays nothing: the mastering chain's bench is flat (3.94 / 3.90 %RT on the heavy case, 4.03 / 4.02 with
every glide moving every call).

New suite `felitronics_saturation_glide_tests`: the snap in every order, a timeline of writes (retarget, shape switch,
a bias glide on Asym, a clock-only pause mid-glide) bit-identical under per-sample, 64, 100, ragged-with-empty-calls,
maxBlock-37 and maxBlock-512 cuts at os 1 and 4 (a mutation that indexes the period from the call instead fails every
row), landing bit-identical to the target stage past the round trip, no click, repeated writes free, no allocation.

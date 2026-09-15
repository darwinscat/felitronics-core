<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam · neural · core — `reset()` restarts the stream, and the price of it is published

`felitronics::nam::NamStage::reset()` was EMPTY, and the comment above it recorded what that cost with a
number: a dense 2001-tap capture that had played a tone answered **digital silence with 0.224604502320**,
and through `reset()` the same. The house verb means a stream RESTART — `eq::EqBand` separated it from
`clearAudioState()` for exactly this reason — so a stale network window speaking into the first receptive
field of the next stream was the defect and not the design. P24 closed the half that is a falling edge (a
lane the host stops handing over is fed the silence it is receiving); this is the other half, the lane that
is PRESENT and gets the caller's own samples.

- **What it now promises is INDEPENDENCE:** nothing the caller fed before the restart can be heard after
  it. Every lane that carried audio is fed the digital silence it still owes, in the call, in full, until
  its state is provably the state of a lane that was silent all along. That promise is **exact** — two
  stages fed different audio before the restart answer the next programme with the same bits, on real
  captures and synthetic ones, at every rate and both widths.
  What it does **not** promise is silence out: a capture answers digital zero with whatever its own biases
  make of it, fresh and restarted alike — measured on NAM's shipped examples at 48 kHz, **0.001195220510**
  for a real Standard and **9.266554832458** for the A2-max feature set. Exact zero is a property of the
  bias-free fixtures, which is what lets them witness the promise.
  Nor, in general, bit-identity with a stage prepared a moment ago: NAM's answer depends on how the stream
  is **cut into calls**, and a restart's chunking is its own. Measured against a stage prepared a moment
  ago: exactly 0 for every fixture in the suite and for a real `slimmable_wavenet`, and **1.037e-06** for a
  real Standard at blocks 64…512, where the restart's last chunk is short — with independence still
  exactly 0 for that same capture, so it is the arithmetic and not the state.
- **The rate-matchers are re-primed too**, because a restart re-anchors the audio-time clocks the way
  `eq::EqBand::reset()` re-anchors its `StateGrid`. Leave the two `core::StreamResampler` legs where the
  previous stream left them and the next programme runs at its sub-sample phase: **1.039e-06 over 5091 of
  5120 samples at 44.1 kHz**, with everything else fixed.
- **NEW API, `core::StreamResampler::clearAudioState()`** — the leg's state alone (zero history, `len`,
  `pos`), keeping the rates, the capacity and the 513 × 64 coefficients. `reset (rates, capacity)` is that
  class's `prepare()`: it reassigns both vectors, `shrink_to_fit()`s on the identity path and re-derives
  every coefficient through a windowed sinc with a Bessel evaluation per tap — **1.77 ms** for one lane's
  two filtering legs (the price is the kernel, not the memory: at an unchanged capacity the allocations are
  reused), and deliberately not `noexcept`. A live stream cannot restart through that; `clearAudioState()`
  is 10 ns. It cannot
  resurrect a refused configuration either: `reset()` leaves `len = 0` when its second allocation throws,
  and a restart that took its length from the buffer's size alone would put that object back to work
  through a coefficient table that is not there, so the class now carries an explicit validity bit.
- **NEW API, `NamStage::clearedSamples()`** — the restart's own odometer, kept apart from
  `drainedSamples()`. The audio cannot witness "the full length for a lane that was playing, the remainder
  for one mid-drain, and NOTHING for a lane that never played": past the debt the output is zero either
  way. A mutation that spends one sample less than the debt survives every audio gate in the suite and is
  caught only here.
- **`NamStage::reset()` is now `noexcept`, and it is an AUDIO-THREAD call whose cost is not the block's:**
  a whole drain length of inference per dirty lane — the field, plus the partitioned-FFT ring any `Linear`
  is charged, plus each rate-matcher leg's tap window, so it is bigger than `prewarmSamples()` and that
  getter is not an estimate of it. On an M-series core, per lane, a real Standard WaveNet is
  **3.77 ms at a 64-sample block — 282 % of that callback** — 3.46 at 256, 3.43 at 512; a real LSTM 1.3 ms;
  a dense 2001-tap `Linear` 0.13 ms. It allocates, locks and throws exactly where `process()` does (nowhere
  for `Linear`/WaveNet; upstream's per-sample Eigen temporaries for LSTM/ConvNet), and it is IDEMPOTENT —
  the debt is re-armed only by audio actually being fed, so a second restart with nothing in between costs
  nothing and a mono host pays for one lane (for a RECURRENT capture it is deliberately not idempotent —
  see below). It also grows faster than linearly as the block shrinks, because NAM's per-call overhead is
  paid `debt / maxBlock` times: a real Standard is 3.61 ms per lane at block 256 and **19.47 ms at block
  1**. There is no cheaper exact mechanism to substitute: NAM's own `Reset` with the prewarm off zeroes
  the Conv1D rings in 0.014 ms and still misses the prepared state by **4089 samples, worst 0.324**, and on
  a `Linear` with the FFT engine it allocates 46 times.
  The allocation carve-out is `process()`'s — the same code path, inherited and not added — and it is wider
  than the header's architecture names suggested: `wavenet_a2_max.nam`, a WaveNet in NAM's own example set,
  allocates **4 times per sample** in `process()` (1024 for one stereo 256-block, measured) and therefore
  inside a restart too. The LSTM/ConvNet row of that carve-out went the other way: at this pin a stereo
  LSTM restart of 48 000 samples allocated **nothing** through either gate, the house operator-new counter
  or Eigen's own.
- **A restart that arrives while the backend is UNPREPARED is parked, not dropped.** `prepare()` writes the
  new `maxBlock` before it can refuse, so an unprepared backend can carry a block of a billion beside a
  256-sample scratch — a restart that touched it is a heap-buffer-overflow, which is why it touches
  nothing there, ledgers included. Re-arming the debt at the next prepare does not cover the lane that is
  PLAYING (its debt is overwritten on every chunk it is fed), so without the parked intent "play, a
  refused prepare, `reset()`, a prepare that succeeds, play" hands back the old stream: **242 samples** of
  a delay(514) capture, measured. The request is honoured at the end of the prepare that can honour it.
- **A recurrent capture stays the NAMED exception, in the mechanism and not only in the comment.** An LSTM
  lane that has already spent its drain reads a debt of zero and is still not empty — the repository's
  slow-cell fixture leaves 0.419413 there — so a recurrent lane that ever played is charged the whole
  half-second heuristic again at every restart (which is what NAM's own `Reset` does) and is never marked
  clean. What that leaves is measured, not promised away: 300 samples differing from a fresh instance,
  worst 1.49e-07, on a real capture.
- **Scope, by byte comparison against v0.33.0:** 6 082 560 samples over six captures × eight rates × three
  block sizes, with loads, clears, mid-stream re-prepares, refused calls, zero-length calls and width
  changes — **byte-identical** where `reset()` is not called. Where it is, 8.394 % of those samples move, by
  up to 0.586079 full scale, and that movement is the fix. On v0.33.0 the same run with `reset()` called
  three times per cell is byte-identical to the run without it, which is the defect stated as a measurement.
- **What the restart does NOT reach, each named with its number** rather than promised away: a recurrent
  cell; NAM's partitioned-FFT clock (swept over nine block sizes x eight rates it peaks at **1.788139e-07**
  against a stage prepared a moment ago, and is EXACTLY ZERO against one clocked to the same point — the
  engine's own arithmetic, not state this stage kept); and a capture
  whose conditioner is a model of its own (`config.condition_dsp`), whose memory neither NAM nor
  `detail::receptiveFieldFromConfig` counts. That last one is a hole in the LEDGER and not in this verb:
  the identical number comes back through the untouched drain — **0.905147969723** after a full drain
  against 0.905148267746 after a restart — and a test now pins both halves of it.
- **Not fixed here, registered:** the conditioner ledger above; `prepare()`, which has the same stale
  window (it is where the 0.224604502320 was first measured); and `rigplayer::RigPlayer`, which has no
  restart verb at all, so a consumer reaching this stage through the player cannot yet call the fix. The
  last two are what a consumer actually hits.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam · rigplayer · convolution — a prepared stage holds no audio, and a player can ask for a restart

P47 made `felitronics::nam::NamStage::reset()` an exact stream restart, and no product could call it.
Two doors stood in the way, and each was half the same defect. This opens both.

- **`NamStage::prepare()` now performs that restart, ALWAYS.** It carried the identical stale window —
  `::nam::DSP::Reset` calls `SetMaxBufferSize` and then a prewarm that is ZERO samples for a `Linear`,
  so `Buffer`'s per-channel input window survived the call. Measured over eight host rates x two block
  sizes x three capture shapes x {re-prepare at the same rate, re-prepare at a different one}, both
  lanes PRESENT: **72 of 96 cells answered digital silence with something, worst 0.567861497402**, and
  0 of 96 do now. A prepared stage and a just-constructed one name one state, so the two verbs of the
  class say one thing.
  There is deliberately **no predicate on what changed**: a re-prepare at the SAME rate and block is the
  common case — a host's buffer-size slider moves more often than its rate one, and a driver stops the
  stream for either — and it is where the leak was loudest (0.468718945980 at 48 kHz, against
  0.469410002232 across a rate change). A "fix it only when something moved" predicate leaves 8 of those
  96 cells leaking, which the mutation stand shows as red.
  The mechanism is P24's ledger and P47's drain, unchanged: `configureRates` charges every lane that has
  EVER been fed a full debt at the new rates, and the tail of `prepare()` spends it. So a first prepare
  after a load costs nothing, and a model change costs nothing either: it prepares a never-fed backend in
  `prepareModel()`, and again in `install()` when the host's numbers moved between the two halves. The `restartOwed_` bit P47 needed
  is **deleted**: once every successful prepare restarts, a parked request has no second question to
  answer, and the sequence it was written for is closed by the stronger rule.
- **THE PRICE, published rather than hidden.** For an architecture whose own `Reset` already prewarms
  (every WaveNet), the drain is a SECOND pass over the field and roughly doubles the call: on an
  M-series core, a stereo real Standard WaveNet `prepare()` on a dirty stage goes **6.5 ms → 13.1 ms**
  at 48 kHz and 7.6 → 14.2 at a 64-sample block; a real LSTM 2.5 → 5.8. That is the message thread, with
  no callback to miss. A first `prepare()` after a load does not move (7.16 → 7.25 ms), which is what
  `everFed_` buys and what the suite gates. Skipping the drain where NAM's own prewarm provably covers
  the ledger is a real optimisation and is registered as one rather than taken here — it is a predicate.
  One consequence is stated because it is a number that moved: against a stage prepared a moment ago, a
  real Standard WaveNet used to be **exactly 0** and is now **1.1e-06**, because the extra drain re-chunks
  a network NAM had already flushed and NAM's answer depends on how the stream is cut into calls. That is
  the same residue `reset()` publishes. INDEPENDENCE — the promise — stays exactly 0 on every capture.
- **`rigplayer::RigPlayer::reset()` — the verb a product can actually call.** The class had none, and
  orbit-amp reaches a `NamStage` only through it (`releaseResources()` is empty there and no host reset
  is overridden). It restarts both model slots (including one the blend law has put to SLEEP), the three
  convolvers **bypassed or not** — a bypassed one is skipped, so its history freezes and is replayed when
  the curve comes back — the dry path's alignment ring, the per-slot whole-sample alignment tails, the
  band filters, the audio-time grid and the scratch; and it SNAPS the gain ramps, as `eq::EqBand::reset()`
  does and for the reason that file records with a number. `prepare()` now calls the same two private
  bodies, so the two verbs cannot drift.
  It deliberately does **not** touch the blend law's state — a restart is not a device change. Re-arming
  the law's warm-up ledger would not deliver its invariant anyway: with `fed = 0` on both slots the law
  ramps its gain down over four blocks, so a by-hypothesis wrong-sounding network is audible regardless,
  and what it buys is a hole — 18 blocks, 192 ms, at every restart. Re-arming only the sounding slot is
  worse: the law rails the goal to the neighbour and plays a full spurious crossfade to the other capture
  and back. The price of the verb is FOUR networks, not one: two stages, each up to two lanes.
- **`convolution::MatrixConvolverNupc::clearAudioState()` (new, additive).** `reset()` there zeroes the
  history AND cancels a swap in flight (`xfadePos_ = 0; state_ = 0`, keeping `cur_`), so a filter
  published a block ago and still crossfading in is dropped — and `CabConvolver::pendingRetry_` is
  already false after a successful publish, so nothing ever re-stages it: the knob move is lost until the
  next knob move. `clearAudioState()` is the history alone, touching only buffers `process()` writes, so
  a composite can restart on the audio thread without losing a filter or racing the message thread.
  `reset()` is now that plus the two cancel lines. `CabConvolver` forwards it.

Green: macOS/clang `ctest` 131/131 with `-DFELITRONICS_WITH_NAM=ON`, deb/gcc-14.2 131/131.

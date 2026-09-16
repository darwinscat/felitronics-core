<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### convolution · lineareq — a restart keeps the operator the caller published, in all three swap-safe convolvers

A swap-safe convolver publishes an operator from the message thread and the audio thread adopts it at the
END of its crossfade, so between those two moments the live slot still names the PREVIOUS operator. All
three `reset()` bodies wiped the publication flag and kept that slot: the operator the caller had already
replaced came back, and nothing was left to re-stage the new one — `setIr()` had already returned true, and
a consumer's retry flag is clear after a successful publish. On those bodies each class's new independence
gate reads the OLD operator in every cell {Pending, Crossfading} x {mono, stereo} x both slot parities
(a worst sample **7.494e-01** `ConvolutionEngine`, **1.017e+00** `MatrixConvolver`, **1.927e+00**
`MatrixConvolverNupc`, from a restart made after the fade settled). The promise those bodies carried was
already the right one — "flush the tail, not revert the EQ" — and reverting to a superseded operator is
exactly reverting the EQ. Law 11e in `docs/DSP-ARCHITECTURE.md` states the rule.

- **Through the product classes** — `lineareq::LinearPhaseEq` and `NaturalPhaseEq` at 48 kHz, a +12 dB
  bell at 1 kHz published over a settled flat curve, the host restarting before a block picked it up — the
  gain at 1 kHz read **+0.00 dB** at 1, 2 and 4 channels (above two channels each channel has its own
  convolver, staged then published) and now reads +12 dB (gated to 0.2 dB). When the curve was the FIRST
  one ever published, the EQ answered **exactly zero** (the fixture's -600 dB floor) until the next band
  move, because the live slot still held the zero taps `prepare()` leaves; it now plays the bell. A restart
  with the move in flight now answers the next programme with the same bits as a restart after it settled
  — a level alone cannot tell `reset()` from an EQ that forwards to `clearAudioState()`, which reaches the
  same curve 20 ms later. `LinearPhaseEq` is also gated through an M/S move (the MSDiag operator, different
  L and R), by level and by independence. The shipped products did not reach this — each re-publishes
  after a re-prepare, and orbit-amp drains silence rather than calling `reset()`, naming this defect as one
  of its two reasons (`src/core/CabinetIr.h:220-222`; the other, the concurrency contract, stays) — so it
  was latent in the products and live in the API.
- **`reset()` now ENDS a swap in flight in favour of the new operator** in
  `convolution::{ConvolutionEngine, MatrixConvolver, MatrixConvolverNupc}`: `cur_` moves to the published
  slot and the state returns to Idle, so `isBusy()` is false at once and the consumer may publish again.
  An operator merely STAGED and not yet published (`stageOperator()` without `publishStaged()`) is
  untouched — it has not been accepted, and the later publish still finds it. The rule is scoped to the
  restart; `prepare()` still discards everything, by contract.
- **WHAT MAKES ADOPTION RIGHT IS LAW 11a INDEPENDENCE, NOT CLICK-FREEDOM.** A half-finished fade is a
  dependency on what came before: two convolvers holding the same published operator, one mid-fade and one
  settled, would answer the next programme differently for up to the length of the fade. Through
  `CabConvolver`, a restart one block into a 50 ms fade and one after it settled now differ in **no** sample.
  **The price is at the seam and it is published** (the engine's suite prints every number here): flushing
  the history is itself a cut — the previous stream's tail stops mid-decay — so the first sample after
  `reset()` steps whichever slot is live, **1.8203e-01** on DC 0.5 into a settled 700-tap operator.
  Adoption puts the new head tap where the old one was at that sample, so it moves the step by at most the
  head-tap difference times the input, in EITHER direction: **7.3203e-01** for a pair whose head tap flips
  +0.70 -> -0.40, **1.3933e-01** for a +1 dB broadband move and **2.2018e-01** for a -1 dB one, and exactly
  the flush for a change that leaves the head tap alone. That difference IS the change the caller asked for.
- **`clearAudioState()` is now on all three** (it was on `MatrixConvolverNupc` alone). It is the history
  and nothing else: it decides nothing about the operator and leaves a crossfade running from where it was
  — which is exactly why it does NOT give independence mid-fade: through `CabConvolver`, a clear one block
  into a 50 ms fade and a clear after it settled differ in exactly the fade's remainder, **2143 of 5120**
  samples a channel (worst 4.186e-01), and the engine's own suite shows 383 of 1600 between two clears
  inside one fade. `rigplayer::RigPlayer::reset()` calls this verb, a choice made when `reset()` still
  dropped the filter; switching it is registered as its own task, not taken here. `reset()` is built from
  `clearAudioState()` in every class, so the two verbs cannot drift.
- **`ConvolutionEngine`'s per-channel `buildIr()` no longer zeroes the staged slot's cached tail.** That
  was a MESSAGE-thread write into a history buffer the audio thread also writes — the same `0.0f`, so
  nothing could be lost, but a data race all the same (put back, ThreadSanitizer reports it), and the one
  thing that made `clearAudioState()`'s promise untrue for this class; its two siblings never write a tail
  off the audio thread. It is redundant: `primeTail()` overwrites that buffer whole at fade start for every
  channel being processed, and a channel NOT being processed was zeroed when it dropped out and is cold,
  for which zero is the right value. Byte-identical output across the change, measured against the tree
  before it over 24 cells of {no reset, reset at Idle} x width {1, 2, 2->1->2, 2->0->2} x block size.
- **The memory order is the one the new write needs.** `reset()` writes `cur_`, which the loader reads
  after acquiring `state_`, so the Idle store is a RELEASE store; and it is made only when a swap was
  actually in flight, because an unconditional `store(0)` would wipe a publication that landed between the
  load and the store. A new real-thread gate, `ConvolutionRestartRaceTests` (outside the wasm-audio tier,
  like core's RtStreams suite), loses **none** of its publications, where the bodies before lose **all** of
  them and the fixed bodies with the store made unconditional lose **12–95 %** per route in every run
  measured (macOS/arm64, Debian/x86-64 with gcc, Windows/x64 with MSVC). With that and the tail write gone,
  `reset()` no longer races a single-producer loader at all: ThreadSanitizer, a loader publishing in a loop
  beside process / reset / clearAudioState on all three classes, reports **0** races against **54–62** on
  the bodies before; the release store put back to relaxed draws **63–71**, and the tail write put back
  **2** (Apple clang needs `-fno-builtin` for the last — it does not instrument a `std::fill` otherwise,
  and a first run here was blind to exactly that; gcc 14 sees it as is). The documented "must not run
  concurrently" contract is kept until a sanitizer row carries that.
- Gates, as a property and not a point: `{Pending, Crossfading} x {mono, stereo} x both slot parities`
  and, for the matrix siblings, every topology with four distinct banks and different L/R inputs, with
  **law 11a independence** as the oracle (two convolvers fed different audio of different lengths, one
  restarted with the publication in flight and one after it settled, answer the next programme — longer
  than the whole FDL span — bit for bit), plus the published operator certified from OUTSIDE the class
  against a reference convolution; `clearAudioState()` mid-fade (still busy, ends on the sample it always
  would have, history gone, lands on the new operator, and not independent of where it was called); for
  the engine, the seam bound and no allocation in either verb; the product classes above; `CabConvolver`'s
  two verbs side by side, the clear by exactly the fade it left running; a fuzz of random publish /
  process / reset / clearAudioState interleavings on all three classes against a settled shadow of the
  same class (with a random width and a reference convolution for the engine), red on the first seed of
  every class before the fix; and the real-thread race gate above. A mutation stand over **17** mutants and
  9 test binaries (the engine, matrix, Nupc, gate-state, `CabConvolver`, fuzz, race, `LinearPhaseEq` and
  `NaturalPhaseEq` suites) is red on every one — including "adopt but leave the state in flight" (the
  operator rolls back one fade later), "`cur_ = 0` instead of `1 - cur_`" (caught only at the second slot
  parity), "leave the FDL", "leave the frame", "forget one slot's cached tail", "the adopted slot's
  topology slips", "the clear restarts the fade", "an EQ forwards `reset()` to `clearAudioState()`", "the
  Idle store made unconditional" (caught by the race gate alone) and each class's pre-P88 body. Four more
  are green on the suites by design: a dropped `ranNc_ = 0` and a dropped `xfadePos_ = 0` are truly
  equivalent; the release store put back to relaxed and the tail write put back are races that only
  ThreadSanitizer sees (above).

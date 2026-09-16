<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### rigplayer · nam — a prepare restates what it counts, and a drained lane is not billed twice

`RigPlayer::prepare()` rebuilt everything the host rate DESIGNS — filters, rings, both stages, the
cold threshold kept in seconds — and nothing the host rate COUNTS. This closes that half, and one
double charge in `NamStage`.

- **A woken slot now warms for the field the NEW rate owes.** The blend law's warm-up debt was latched
  from `warmFor()` when a model landed and read for ever after, including by the wake of a sleeping
  slot. Measured through the player on a 6x6 host-rate grid, two block sizes and three capture shapes:
  the warm-up after a rate change did not depend on the new rate at all. **48 → 96 kHz warmed a slot
  for 2048 host samples where 4096 are owed — half the field; 48 → 192 kHz a quarter; 192 → 48 kHz held
  a slot silent four times too long.** Of 180 cells, 85 under-warmed, 85 over-warmed and 10 read right
  — eight of those 44.1 ↔ 48 kHz, so the pair a fixture reaches for first could not see it. **180 of 180
  now match a player prepared at the new rate from the start.**
- **The ledger is restated by the rule each count's algebra allows**, in a new law verb,
  `nam::blendRestated()` — not `blendLanded()`, which would also clear a load in flight, a sleep and a
  refusal. The debt is RECOMPUTED: it is not homogeneous in the rate (the rate-matcher's latency is 0 at
  48 kHz, 61 at 44.1, 96 at 96, and the block term does not scale), so rescaling it is a second copy of
  `warmFor()` with a different answer. The progress is mapped by its PREDICATE: an audible slot stays
  audible, a warming slot starts over. Scaling it by the rate ratio is the trap — 96 → 48 kHz turns a
  just-audible slot inaudible for any block over 96 samples and plays a spurious full crossfade. The rest
  count IS rescaled: it is pure elapsed time, so the ratio is exact, and is 1 where nothing moved.
- **A restart inside a warm-up re-arms it, at any rate — `reset()` included.** Both verbs flush every
  network, and the law went on crediting a warming slot the field it heard before the flush. Restarting
  1, 4 and 8 blocks into a 2001-tap field left 14, 11 and 8 blocks of warm-up against the 15 a restart
  at the landing costs. A warming slot is at weight zero, so re-arming it is silent; an AUDIBLE slot is
  still not re-armed, for the reasons P85 gave.
  *This supersedes the P85 note above that `RigPlayer::reset()` does not touch the blend law's state:
  it does not, except for this.*
- **The rest a slot has served survives a restart.** It was left in the old rate's samples against a
  threshold recomputed in the new: 48 → 96 kHz slept 0.73 s late, 96 → 48 kHz after one block. A first
  draft zeroed it instead — a behaviour change at the UNCHANGED rate, postponing every parked dial's sleep
  by two seconds at each restart — and was caught before it shipped.
- **The per-slot alignment delays are restated too** — they are host samples, written only by a knob, a
  table or a landing — and **a landing delivered before a rate change and taken after it** warms and
  aligns at the new rate. A `prepare()` between `load()` and the next block leaves the previous pack's
  model ids alone: they index a model list that now belongs to the new pack.
- **Slot trims switched OFF stay off across a restart.** The snap ignored `setInputTrims(false)`, so
  a restart ducked the slot to the pack's trim and ramped back over ~43 ms (−4.29 dB on the first block
  of a −6 dB entry). The shared body fixes `prepare()` and `reset()` at once.
- **`NamStage`: a lane whose falling-edge drain ran to the end is not billed again by the next
  `prepare()`.** Only a restart cleared its "may be holding audio" flag, so the next prepare charged a
  whole drain for a clean lane: 4093 samples for `wavenet_a1_standard`, 6347 for `A2`, and the
  re-prepare measured 14.7 ms where 12.7 is owed. A PARTLY drained lane is still charged in full, and a
  recurrent one always is. *This supersedes the P85 note above that `configureRates` charges every lane
  that has EVER been fed.*
- **`rigplayer` has a permanent RT allocation gate** (`RigPlayerRtAllocTests`), on the shared counter
  that sees all eight forms of `new`: `process()` on five shapes, across a turn, through a slot falling
  asleep and draining, a host dropping to mono and to a gap, a call four times the block; `reset()`; the
  per-block read-outs. Its control plants an allocation and must read one more than the same window
  without it.

**Registered, not taken:** skipping `prepare()`'s second warm-up pass for architectures whose `Reset`
already clears their state (P98). With the drain removed outright, every capture NAM ships keeps
independence at exactly 0 — and this tree's own `Buffer`-based fixtures leak on 72 of 96 cells, worst
0.567861. NAM's example set has no such capture in any tree, so a check against real captures alone
would have approved a blanket skip. Found on the way and registered: the alignment ceiling
`kBlendMaxDelay` is 128 samples, so a 64-sample lag at 48 kHz is clamped to half at 192 kHz (P99).

Independence (P47/P85) is unchanged: exactly 0 on every capture.

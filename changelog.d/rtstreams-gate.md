### core — `RtStreams` closes its two races: the audio thread holds a busy bit, the controller's calls refuse instead of writing under it

**BREAKING (API).** `AuditionStream`, `ConvStream` and `RecStream` changed shape on purpose: an audio callback
written for the old types no longer compiles, rather than compiling and losing publishes.

Two races are gone (#273):

- **A publish under a playing block.** `AuditionStream::publish` (and `ConvStream::refill`) could overwrite the
  samples a callback block, already past its flag check, was still reading. The block could also store its old
  position and its end-of-clip flag after the new clip was armed, so the new clip started midway or never played.
  ThreadSanitizer reported it; the CI job carried a suppression for it, now removed.
- **A recorder rewind.** A push whose `recording` check came before `stop()` could land after the next `start()`:
  either at `buf[0]` of the new take, unordered with the previous harvest's read, or by storing the OLD length + 1
  over the rewind, so the new take began with the old take's samples. CI's TSan job caught it intermittently.

**The gate.** `RtGate` is one atomic word per stream: the mode in the low 31 bits and a BUSY bit the audio thread
holds while it is inside a block it was admitted to. The word changes only by read-modify-write. A controller call
takes the mode down and learns, in the same operation, whether a block is open. If one is open, it refuses and
changes nothing; if none is, it mutates and raises the mode again.

- There is no counter and no wait. A stopped device holds no block, so every call lands at once.
- A refused call lands on a retry once the block closes. The block's scope is all the retry waits for, and nothing
  can open another block while the mode is down.
- A Block that never closes keeps BUSY set: every later call is refused, and the stream is silent rather than torn.

**The API.**

- **Audio thread:** `auto b = s.beginBlock();` is one RAII `Block` per stream per callback. It is empty unless the
  stream plays, and it reads and advances the stream only while it is open. `finish()` ends a one-shot clip at
  close. `RecStream::Block::push()` appends.
- **Controller** (one thread, sizing included): `[[nodiscard]]` `tryReserve`, `tryPublish`, `tryRefill`,
  `trySetMode`, `tryStart`, `tryHarvest(std::vector<float>&)`, plus `stop()`. A refusal mutated nothing and left
  the stream stopped; keep the data and retry from a timer.
- **Any thread:** `isPlaying()`, `mode()`, `isRecording()`, `recorded()` and `AuditionStream::loop`.
- **Three behaviour fixes:**
  - **`ConvStream` keeps its clip in FIXED storage with its own length.** `refill` no longer reshapes the vector the
    audio thread reads.
  - **A route switch cannot be taken down by the old clip's end.** Live monitoring used to switch off silently when
    a clip block ended it just after the switch.
  - **`RecStream::tryReserve` discards the take.** The old `reserve` kept a length past a shrunk buffer.
- **Capacities are clamped to `INT_MAX`**, never asserted. A false from a try therefore always means "a block was
  open".

`RtGate` is public, so a product can gate its own swaps of plain memory the same way.

**Proof.**

- `felitronics_core_rtstreams_interleave_tests` runs one controller call at every atomic operation of an audio
  block, and one block at every operation of a controller call, on one thread through a hooked atomic. It runs in
  every tier, wasm included.
- `felitronics_core_rtstreams_tests` forces contention with real threads: the audio thread holds a block open, the
  controller must be refused, then must land. The CI TSan job runs it with nothing suppressed.
- Every protocol mutant turns at least one of the two suites red (see the PR).

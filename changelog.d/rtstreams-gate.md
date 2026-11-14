### core — `RtStreams`: the audio thread holds a busy bit, and the controller's calls refuse instead of racing

**BREAKING (API).** `AuditionStream`, `ConvStream` and `RecStream` changed shape on purpose: an audio callback
written for the old types no longer compiles, rather than compiling and losing publishes.

Two races are gone:

- **A publish under a playing block (#273).** `AuditionStream::publish` (and `ConvStream::refill`) could overwrite the
  samples a callback block, already past its flag check, was still reading. The block could also store its old
  position and its end-of-clip flag after the new clip was armed, so the new clip started midway or never played.
  ThreadSanitizer reported it, and the CI job carried a suppression for it, now removed.
- **A recorder rewind** (CI's TSan job found it while #273 was open). A push whose `recording` check came before
  `stop()` could land after the next `start()`: either at `buf[0]` of the new take, unordered with the previous
  harvest's read, or by storing the OLD length + 1 over the rewind, so the new take began with the old take's
  samples.

**The gate.** `RtGate` is one atomic word per stream: the mode in the low 31 bits and a BUSY bit the audio thread
holds while it is inside a block it was admitted to. The word changes only by read-modify-write. A controller call
takes the mode down and learns, in the same operation, whether a block is open. If one is open, it refuses and
changes nothing but the stop; if none is, it mutates, and raises the mode again if the call starts something
(`tryPublish`, `trySetMode`, `tryStart`; `tryRefill`, `tryReserve` and `tryHarvest` leave it down).

- There is no counter and no wait. A stopped device holds no block, so every call lands at once.
- A refused call lands on a retry once the block closes. The block's scope is all the retry waits for, and nothing
  can open another block while the mode is down.
- A Block that never closes keeps BUSY set: every later call that needs the stream to itself is refused (`stop()`
  and `trySetMode(0)` still take it down), and the stream is silent rather than torn.

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
  - **A `tryReserve` that lands stops the stream and empties it,** and resets its length BEFORE it allocates: a
    reserve whose allocation throws leaves an empty stream, never a length over freed storage. (A refused one only
    stops the stream, like every refusal.) The old
    `RecStream::reserve` kept a length past a shrunk buffer.
- **Capacities are clamped to `INT_MAX`**, never asserted. A false from a try therefore always means "a block was
  open".

`RtGate` is public, so a product can gate its own swaps of plain memory the same way.

**Proof.**

- `felitronics_core_rtstreams_interleave_tests` runs one whole controller call at every atomic operation of an
  audio block, and one whole block at every operation of a publish, a refill and a start, on one thread through a
  hooked atomic. Where exceptions exist it also tests a reserve whose allocation fails. It runs in every tier,
  wasm included.
- `felitronics_core_rtstreams_tests` forces contention with real threads: the audio thread holds a block open, the
  controller must be refused, then must land. It reserves under a live reader, and it forces the ABA case where the
  admitting CAS succeeds on a re-armed word. The CI TSan job runs it with nothing suppressed.
- Each of 27 mutants of the protocol turns at least one suite red. Among them: a skipped wait in every try, a close
  that publishes before its stores, an arm before the copy, every relaxed order (TSan's model kills those; on x86
  the hardware cannot tell them apart), a reset after the allocation. Two weakenings the protocol does not rely on
  pass, as they should: a relaxed load before the admitting CAS, and an acquire-only admitting CAS.

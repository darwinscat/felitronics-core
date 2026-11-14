### delay — `TapeDelay` reads its tape with up to four playback heads

One head on the line meant every repeat was the one before it passed through the loop again, so the repeats
always fell: no setting made the third echo louder than the second. Now up to `kMaxHeads` (4) heads read the one
tape (#272).

- **A head** (`TapeDelayHead`) has a POSITION, a fraction of TIME in (0, 1], a LEVEL (0 … 1, 0 is off) and a PAN
  (−1 … +1). TIME stays the one motor: a head stands at its position times the motor's glided time, so turning TIME
  moves every head, and a head moved on its own glides under TIME's law and bends that repeat, as a tape echo does
  when a head is shifted. No head plays past TIME, so the line keeps its length and memory.
- **The loop stays at 1x TIME.** REPEATS feeds back from the full length, never from a head, so the heads' figure
  repeats once per TIME and decays with REPEATS: heads at 1/3, 2/3 and 1 play a triplet; at 1/4, 2/4, 3/4 and 1,
  four on the beat. Moving a head changes the figure, not the loop's rhythm. OFFSET is unchanged: one time skew
  between the channels.
- **Pan is a balance** on the head's own lines, at constant power with the centre at unity on both sides, so a
  hard-panned head is +3 dB on its side. The law is rational — gL = (1 − p)/√(1 + p²), gR = (1 + p)/√(1 + p²),
  exact at the centre and the ends, the same bits on every row with no transcendental — and within 4.1° of the
  sine law. On dual mono (one guitar on both lines, as the amps and the demo feed it) heads panned apart are a
  ping-pong; on a stereo input a hard-panned head plays its own side only. A mono call has no pan: every head plays
  at its level.
- **Moved live.** Level and pan reach a head as two side gains on counted 20 ms ramps that land exactly; a silent
  head stands where it was sent, so turning one on plays it there instead of sweeping it in; a change of call width
  ramps instead of stepping. All of it is sample-clocked: the same playing gives the same samples whatever the
  host's blocks, live head moves included.
- **API.** `TapeDelayParams` gains `heads` (a `std::array` of `kMaxHeads`), applied by `setParams` after the five
  knobs: the struct is the whole state, and one built without heads carries the default layout. `setHead (k, head)`
  and `setHeads (heads)` set them directly; each field is bounded on its own and an index off the array is ignored.
  `shownHeadMs (k)` publishes where each head reads, for the picture.
- **The default is today's sound, byte for byte**: one head at 1x, level 1, centre. Every existing preset, setter
  call and `setParams` without heads plays as before — held to a frozen copy of the one-head class by a parity
  test that compares the output and every tap as bytes through prepares, setters, every call width, stand-by,
  resets and poison.
- **`tailSeconds()` counts the heads**: stacked heads add and hard pan is +3 dB, so the first echo is taken at the
  sum of each head's level times its louder side. The default layout reports exactly what it did.
- **MIX** still adds the wet with the dry untouched; with several heads the wet is their sum, so at MIX 1 four
  stacked unity heads play four times one.

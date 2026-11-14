### convolution — `MatrixConvolverNupc`'s last stage no longer lands in one call

The last stage of the schedule (blocks of 2048) holds almost every partition of a long IR, and it did all of them
at its boundary: one `process()` call in every 2048 samples did the whole stage's spectral MACs. For a 20 s
true-stereo (Full) IR at 48 kHz that is 1872, and a swap's first call primed every stage of the incoming slot on top
(#270). The worst call outgrew a 64-sample block on an Intel i9 core and, in a swap, a 128-sample one.

Now that stage is **spread** for an operator whose deferred work there — (P − 1) partitions times the banks its
topology reads — reaches `kSpreadFromMacs` (256 MACs a hop).
- **The schedule.** For the hop ahead, every product with j ≥ 1 reads history that is already final. Those products
  are accumulated through the hop on a schedule of sample positions, so the work and its order never depend on how
  the host slices the stream. The boundary adds only j = 0 and inverts.
- **The swap.** A spread operator picked up over a warm history is prepared for one hop while the live one plays
  alone, and its fade starts one hop (2048 samples, 43 ms at 48 kHz) after the pickup. Over a cold history (nothing
  through that stage since `prepare()` or `clearAudioState()`) the fade starts at once.
- **Output.** A spread operator's sum is taken in another order, so its output differs from the unspread one by
  float noise. Below the threshold everything runs as before, to the bit and with an immediate fade: at 44.1 and
  48 kHz that is every cabinet and every lineareq EQ in the tree.
- **New API.**
  - `maxSwapSamples()` (also on `CabConvolver`): the most samples a swap can take, the crossfade plus a hop where
    the geometry can spread.
  - The threshold is the class's second template parameter, `MatrixConvolverNupc<Fft, SpreadFromMacs = 256>`, a
    test seam.
- **Memory.** `prepare()` adds the spread state, 96 KiB at most, only where a cap-length operator could spread.

**cabinet.** `CabinetIr` resting out of the path fed a handed-over IR 50 ms of silence and then stopped. A spread swap
was then still in flight, and the cabinet that came back played the rest of its fade. It now feeds what the
convolver's `maxSwapSamples()` says.

**test_support.** The allocation counter also counts releases (`alloc::frees`). A request count cannot see a buffer
kept that should have been freed; the convolver's mono re-prepare is gated on it now.

**ci.** A ThreadSanitizer job (`FELITRONICS_ENABLE_TSAN`) runs the suites with real threads. Among them is a new one
where a loader publishes growing, shrinking and topology-changing operators beside an audio thread.

**tools.** `fcore_convburst` measures the worst call against the block's budget, out of a swap and through one, with
its 99.9th percentile and how long a swap keeps the convolver busy.

Measured with `fcore_convburst`: worst `process()` call, 20 s Full IR under a 20 s cap at 48 kHz, 50 ms fade, pffft.

| machine | block | idle, before → after | swap, before → after |
|---|---|---|---|
| Apple M5 Pro | 64 | 64 % → 6.3 % | 190 % → 12.7 % |
| Apple M5 Pro | 128 | 33 % → 4.8 % | 105 % → 10.0 % |
| Intel i9-13900H (unpinned) | 64 | 174 % → 10 % | 460 % → 18 % |
| Intel i9-13900H (unpinned) | 128 | 77 % → 12.5 % | 245 % → 18 % |
| Intel i9-13900H (E-core) | 128 | — → 26 % | — → 36 % |

A swap of such an IR keeps the convolver busy 93 ms (the 50 ms fade and the one-hop preparation); the 2 s IR is under
the threshold, swaps in 50.7 ms as before, and its worst call is 19 % of a 64-sample block on the M5, 34 % on the i9.

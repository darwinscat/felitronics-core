### convolution — `MatrixConvolverNupc`'s last stage no longer lands in one call

The last stage of the schedule (blocks of 2048) holds almost every partition of a long IR, and it did all of them
at its boundary: one `process()` call in every 2048 samples did the whole stage's spectral MACs. For a 20 s
true-stereo (Full) IR at 48 kHz that is 1872, and a swap's first call primed every stage of the incoming slot on top
(#270). The worst call outgrew a 64-sample block on an Intel i9 core and, in a swap, a 128-sample one.

Now that stage is **spread** for an operator whose deferred work there — (P − 1) partitions times the banks its
topology reads — reaches `kSpreadFromMacs` (64 MACs a hop).
- **The schedule.** For the hop ahead, every product with j ≥ 1 reads history that is already final. Those products
  are accumulated through the hop on a schedule of sample positions, so the work and its order never depend on how
  the host slices the stream. The boundary adds only j = 0 and inverts.
- **The swap.** A spread operator picked up over a warm history is prepared for one hop while the live one plays
  alone, and its fade starts one hop (2048 samples, 43 ms at 48 kHz) after the pickup. Over a cold history (nothing
  through that stage since `prepare()` or `clearAudioState()`) the fade starts at once.
- **Output.** A spread operator's sum is taken in another order, so its output differs from the unspread one by
  float noise (about 1e-7 of the peak). Below the threshold everything runs as before, to the bit and with an
  immediate fade: an operator under 17 partitions of that stage in Full (34816 taps), 33 in LRDiag or MSDiag
  (67584 taps: a cabinet up to 1.41 s at 48 kHz, 0.70 s at 96 kHz), 65 mono. lineareq's top quality (131073 taps)
  now spreads in stereo, Full and LRDiag / MSDiag alike: a band move lands a hop (43 ms at 48 kHz) later, beside the
  1.37 s of latency that quality carries. Its other qualities, its per-channel lanes on a bus wider than stereo
  (63 MACs a hop at the top quality) and NaturalPhaseEq stay on the old path.
- **A consumer that waits a fixed fade for a swap to land must ask `maxSwapSamples()` instead.** A cabinet IR longer
  than 67584 taps now starts its fade a hop after the pickup. Core's `CabinetIr` asks; a copy that still waits 50 ms
  (orbit-amp's own) comes back mid-fade after a handover made while it rested.
- **Why 64.** Below it a stage's deferred work is too small to be worth a hop of delay. Just under it the old burst
  stays, about 3 × 63 MACs in an unspread swap's first call — per convolver the worst call left in the tree, and a
  small one (measured below). A composite that publishes several convolvers back to back (lineareq's lanes on a bus
  wider than stereo) pays it once per lane in the same callback, as before.
- **New API.**
  - `maxSwapSamples (numBanks = 4)` (and `CabConvolver::maxSwapSamples()`, which asks for its two banks): the most
    samples a swap can take, the crossfade plus a hop where an operator of that many banks could spread at this cap.
  - The threshold is the class's second template parameter, `MatrixConvolverNupc<Fft, SpreadFromMacs = 64>`, a
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

Measured with `fcore_convburst`: the worst `process()` call as a share of the block's budget (its 99.9th percentile in
parentheses), Full IR under a 20 s cap at 48 kHz, 50 ms fade, pffft. "Before" is v0.64.0.

| machine | IR | block | idle, before → after | swap, before → after |
|---|---|---|---|---|
| Apple M5 Pro | 20 s | 64 | 64 % → 6.3 % (4.6 %) | 190 % → 12.7 % (10.1 %) |
| Apple M5 Pro | 20 s | 128 | 33 % → 4.8 % (4.2 %) | 105 % → 10.0 % (9.0 %) |
| Intel i9-13900H, unpinned | 20 s | 64 | 174 % → 33 % (9 %) | 460 % → 41 % (19 %) |
| Intel i9-13900H, unpinned | 20 s | 128 | 77 % → 19 % (8 %) | 245 % → 23 % (21 %) |
| Intel i9-13900H, E-core | 20 s | 128 | — → 26 % (16 %) | — → 37 % (32 %) |
| Intel i9-13900H, E-core | 2 s | 64 | — → 10 % (8 %) | 69 % → 16 % (14 %) |
| Intel i9-13900H, E-core | 34816 taps, just under | 64 | — → 17 % (13 %) | — → 31 % (30 %) |
| Intel i9-13900H, E-core | 36864 taps, just over | 64 | — → 10 % (8 %) | — → 16 % (13 %) |

A spread operator's swap keeps the convolver busy 93 ms (the 50 ms fade and the one-hop preparation); one under the
threshold, 50.7 ms as before. The i9's isolated maxima of about 0.4 ms against a far lower 99.9th percentile recur on
idle rows too and are the shared machine, not the code: the work per call is bounded by the test suite's counts.

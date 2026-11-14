### convolution — `MatrixConvolverNupc`'s memory follows the loaded IR

`prepare()` used to size every bank spectrum to the cap, for four banks in each of the two slots whatever the
topology: about 80 bytes per sample of cap in stereo, most of it never read by a short IR, and half of it never read
by LRDiag or MSDiag. The bank spectra now follow the IR: `stageOperator()` (and `setOperator()` / `setIr()` through
it) grows each bank, stage by stage, to the partitions the IR it is handed reaches, and only the banks the topology
reads (mono 1, LRDiag and MSDiag 2, Full 4). They grow only — a shorter IR reuses what a longer one left — so each
slot holds its own high-water mark until the next `prepare()`, which releases it. A re-prepare to a smaller cap (a
lower host rate) or to mono now gives the history's memory back as well; `assign()` used to keep the old capacity,
and a mono re-prepare kept the right channel's whole history.

The shared history (the FDL) stays sized to the cap, on purpose: both slots read that one history, so a new operator
answers the last cap of input at once. The swap stays click-free, and for every sequence of accepted calls the
output is what it was, bit for bit.

What changes for a caller: `setOperator()` / `stageOperator()` allocate when an IR reaches further than that slot
has held before (on the message thread, where the class already documented that they allocate). `process()` still
allocates nothing. Exhaustion during that growth throws out of the call (law 11d); a native caller that catches it
outside any `noexcept` frame keeps an instance that plays on and builds the next swap correctly, and an operator it
had staged before stays intact. Through a `noexcept` caller, or on the wasm tier, exhaustion is fatal, as the law
says. A swap refused for a null bank is now refused before it touches the inactive slot; it used to relabel the slot
and rebuild the banks ahead of the null one first, so publishing after such a refusal now publishes what was staged
before it.

New, `reserve (topology, len)`: grows both slots' banks now, so that no later swap up to that length allocates. It has
`prepare()`'s thread contract. `CabConvolver` (the amps, OrbitCab) gets the new behaviour with no change of its own.
lineareq's EQs call `reserve()` from their `prepare()` for every bank either slot can hold: their IR is always the
cap's length and any band set may need the Full matrix, and their `setBands()` promises to allocate nothing (it is
`noexcept`). Their memory and their allocations are exactly what they were.

The class's COST NOTE claimed that CPU is set by the cap, and that was overstated. The input writes, the head and each
stage's FFTs follow the schedule: with a 128-sample head that is five stages for any cap above 2048 samples. The
spectral MACs, the bulk of a long IR's work, follow the loaded IR (`activeParts`). The note now says so, beside a
MEMORY NOTE stating the same split for memory.

Measured (`modules/convolution/tests/MatrixConvolverNupcMemoryTests.cpp`, scalar and pffft backends):
- **Bit-exact against the class it replaces.** Compared with a frozen copy of the previous class, every output
  sample is identical, fades included. That holds through sequences short, long, short, long, past the cap, cleared
  and back, for mono, LRDiag, MSDiag, Full and all of them mixed in one instance, at uneven blocks and at blocks of 1,
  64 and 4096. It also holds with a swap tried mid-fade, `reset()` and `clearAudioState()` mid-fade, and an operator
  staged and published late.
- **The full IR after growth.** A cap-length IR loaded after a short one is the direct convolution of the whole
  input, within 5e-7 of a 1.7 peak, and its last tap answers an impulse.
- **Bytes against the schedule.** What each swap requests equals an oracle rebuilt from the schedule's rule; a
  shorter or equal IR requests nothing; `prepare()` requests exactly the old bytes minus the bank spectra.
- **Exhaustion.** Each of the 16 requests a growth makes was failed in turn, twice over. Each throws. The instance
  then plays on and, after the next swap, bit-exact against one that never failed; an operator staged before the
  failed swap publishes intact.
- **`reserve()`.** It requests both slots' banks, the oracle's size, and six swaps of every topology after it
  request nothing. lineareq's `setBands()` allocates nothing across MSDiag, LRDiag and Full at 1, 2 and 4 channels.
- **Nothing in `process()`.** No allocation in `process()` anywhere in the suite.

Bytes per stereo instance, both slots loaded, at 48 kHz: what the class requests from `operator new`, scalar backend
(pffft's own plan allocations go through `malloc` and are not counted). The test prints 44.1, 48, 96 and 192 kHz, for
both backends.

| case | before | after |
|---|---|---|
| amp, cap 5 s, LRDiag, 1.28 s IR | 19.7 MB | 6.2 MB |
| cab, cap 4 s, LRDiag, 1.28 s IR | 15.7 MB | 5.4 MB |
| reverb, cap 20 s, Full, 2 s IR | 77.2 MB | 21.9 MB |
| reverb, cap 20 s, Full, 4 s IR | 77.2 MB | 28.0 MB |
| reverb, cap 20 s, Full, 20 s IR | 77.2 MB | 77.2 MB |

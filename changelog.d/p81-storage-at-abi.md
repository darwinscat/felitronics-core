### tools · analysis — the price of a measurement, asked before it is paid

`fc_probe_report_storage_bytes`, `_bursts_`, `_hum_`, `_forensics_` and `_lowend_storage_bytes` publish,
through the wasm ABI, the number each offline analyzer already computes for itself before it allocates a
byte (law 11d). Until now nothing carried it out of the module: JavaScript called `_run`, which went
straight into an allocating `prepare()`, and there was no way to ask the price first. The prices are not
small — `HumDetector` at 768 kHz and 16 channels asks for 352 688 184 bytes, the five together for
377 423 824, against a linear memory that stops at 2 GiB and also has to hold the caller's decoded input.

Each query takes the GEOMETRY only — `(channels, sampleRate)` — because the whole use is to ask before the
input buffer exists; `_run` keeps its own buffer checks, and a price above zero is not a promise about a
particular call. A refused geometry quotes a canonical `+0.0`. The answer is a `double`, not the
`std::uint64_t` the budget is accounted in: an i64 does cross this boundary on the pinned toolchain, as a
BigInt, and a BigInt is the wrong shape for a number a page has to add to its own input size and put in a
report — `bigint + number` throws and `JSON.stringify` refuses it, so every use site would need a
conversion first.

What the module does today when the allocation fails was measured rather than assumed, and the assumption
was wrong in both directions. It is worse than a bad return value: `-fno-exceptions` turns the failed
`operator new` into `abort()`, so the C entry point never returns at all — the caller reading
`_fc_probe_hum_run(...) === 1` is handed a thrown `WebAssembly.RuntimeError` instead of a status. It is
also not the death of the page: that error is catchable, and after catching it the module went on
answering, with the aborted mode's getters reading 0 as the `have*` discipline promises. No ceiling on
geometry was added — what a page can afford is the page's number, not this file's.

Also in this change, found while building the gate for it:

- **The export whitelist in `tools/wasm/build.sh` had four ways past it**, none of which its guard could
  see: it recognised only a return type of `int|double|std::uint32_t`, it matched the declaration at
  column zero, and its guard was a hand-written floor (`-ge 39` against 86 actual entry points). So an
  indented declaration was invisible to every scanner at once while `EMSCRIPTEN_KEEPALIVE` exported it
  anyway; two declarations on one line dropped the second from the list with the counts still agreeing;
  a one-line body calling another `fc_*` function exported the wrong symbol; and a comma declarator
  (`FC_EXPORT int fc_a (void), fc_b (void);`) is one token, one line and one name while the second is
  simply absent. All four were reproduced against the gate before they were fixed. Both ABIs' lists are now read as the identifier
  before the first `(`, one declaration AND one declarator per line are enforced, and the count is
  checked for EQUALITY against the declarations; a second gate refuses a return type the boundary does not carry rather than
  dropping it. `fc_master`'s copy had the type defect too, plus a name pattern without digits that would
  have TRUNCATED `fc_render_v2` into the list rather than omitting it.
- **Two rationales in the tree were false** and are corrected with the measurement that settled them:
  `-sWASM_BIGINT` is ON by default in emscripten 6.0.9, and `EMSCRIPTEN_KEEPALIVE` does keep a function
  that is missing from `-sEXPORTED_FUNCTIONS`.
- **`planarSpan()`'s channel bound was pinned by nothing in the repository.** Deleting it left all 123
  tests green, because every other entry point refuses the width a second time in `prepare()` — every one
  but `fc_probe_needle`, which has no `prepare()`. `felitronics_abi_tests` now pins it there.
- The five modes' default parameters are one `constexpr` constant each, read by the run, by the price and
  by the result getters that used to construct their own.

The gate: `felitronics_analysis_abi_tests` goes from 389 checks to 575 — the allocation a first run
actually asks for against the published demand, measured at the WIDEST geometry the ABI has and equal to
the byte; the price against the core's own `storageFor()` at every width and rate, each row asserted
positive so it cannot pass as `0 == 0`; a canonical `+0.0` on 75 refused geometries; a grid of five modes
by twenty-six rates by nine widths by three programme lengths, where the price must be positive exactly
where `_run` is accepted and every refused row must leave its getters silent; and the asymmetries that
are deliberate, written down so nobody "fixes" them. `tools/wasm/storage-probe.mjs` makes the assertions
only the wasm tier can make — starting with the one no native test can, that all five names reached the
artifact — and prints the demand table that `felitronics_analysis_abi_tests --storage-table` prints
natively: 425 rows, byte-identical across the two tiers, on the release and the checked module alike.

Twenty-six mutants were run against it and twenty-three died. **The crew's testing round found the hole
the first twenty missed**: with `_run` ignoring what `prepare()` returned, the whole suite stayed green
while an empty programme at a refused geometry was ACCEPTED and its getters served the previous
programme's numbers — `_run` skips `process()` when there are no frames, so an unprepared analyzer still
reached `finish()`. Both grids had a single non-zero programme length. They now walk an empty one, a
short one and a full one, in both tiers, and the mutant of that line dies in four of the five modes —
in lowend it is equivalent, because `planarSpan` refuses an empty programme before `prepare()` is reached. Two more
survivors from the same round — an allocation that only happens above two channels, and one that only
happens away from 48 kHz — are what moved the allocation oracle to the widest geometry. The three that
survive are equivalent and named where they live: removing `setParams (kReportParams)` from report's run,
where the constant holds the instance's own defaults; hard-coding the width in the lowend query, whose
demand does not depend on it; and masking the low three bits off every price, since every demand in the
accepted domain is a multiple of eight. Six mutants of the build gate were killed too — the indented
declaration, the two on one line, the comma declarator, the empty extraction, the unsupported return
type, and the name taken from a body call.

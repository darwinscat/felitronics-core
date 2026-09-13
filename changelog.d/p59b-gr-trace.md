### `mastering` · `tools` — the gain-reduction trace: WHERE the compressor and the limiter worked (ABI v4)

A loudness search used to report its gain reduction as three numbers for the whole programme — mean, p95, max. The
solution now also carries each stage's TRACE: the delivered programme cut into uniform buckets, and per bucket the
largest and the mean |GR|.

- **`mastering::GainReductionTrace`** in `LoudnessSolution` (`compressorTrace`, `limiterTrace`), built by
  `GainReductionTraceBuilder` in the solver's tap sink — the same tap values, windows and units the statistics use
  (frames for the compressor, frames × `tapOversampleFactor` sub-samples for the limiter), so where the statistics
  are valid the maximum over the buckets IS their `maxDb`, bit for bit. A maximum per bucket, not a point sample.
  `min(1000, frames)` buckets over the programme's frames, boundaries `floor(k·F/B)`, never an empty bucket; per bucket
  a sample count and a non-finite count, so a poisoned or empty stretch cannot read as "idle".
- **It describes the audio in `out`.** Reset and rewritten on every render, including the delivery re-render of a
  search that did not end on its best point; `valid` only when the render ran to its end, the window saw a sample and
  none was non-finite; no buckets for a verdict reached before any render. It lives in the solution, so a later solve
  does not touch it. The price is the solution's size: 2336 → 50 384 B per solution record, pinned with its formula.
- **ABI v4** (the rule's third bump): the new entry point `fc_solution_gr_trace(solution, stage, out, cap, written)` on
  the pattern of `fc_solution_log`, the header-less `fc_gr_trace_bucket` (frozen from v4), the `fc_gr_stage` codes,
  and the traces' bucket counts and validity at the end of `fc_measurement` — one row of the size table (224 B from v4).
  `fc-master-layout.mjs` is at v4; the site's copy of it must follow before a page reads the trace.
- **Verified:** the trace nulls bit for bit against the chain driven by hand and bucketed another way (three lengths,
  including one shorter than 1000 frames); the ABI's traces are the core's bit for bit — on a plain handle with both
  stages live (MasterAbiTests) and on the delivered-rate path (`fcore_master selftest`, whose fixture leaves the
  limiter idle) — and native agrees with wasm through `master-parity.mjs`; the re-render branch, the no-re-render
  branch, an early exit, the refusals and a failed first render each carry the trace their render left; an impulse's
  GR sits in its own bucket and ends where the limiter's hold plus release puts it. Built and run with MSVC as well.
  Mutation stand: 22 mutants, 18 killed; the four survivors are equivalent on every reachable input (the two stages'
  bucket counts and validity are equal by construction, the one reachable `RenderFailed` refuses before any tap, and a
  defensive branch for a stream that goes backwards).

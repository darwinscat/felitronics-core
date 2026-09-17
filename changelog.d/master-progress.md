<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — progress and cancellation of the whole-programme calls; C ABI v5

**`felitronics/mastering/Progress.h`.** `ProgressCallback { bool (*fn) (void* context, const ProgressEvent&); void* context; }`,
taken by new overloads of `TargetLoudnessSolver::solve` / `measureInputLoudnessRange` and `DeliveredMastering::solve` /
`measureInputLoudnessRange`. An event carries the stage (`Convert`, `LoudnessRange`, `SearchPass`, `FinalRender`), the
render's number and its bound, the fraction of the stage — exactly 0 first, exactly 1 last, no two events more than 1 %
of the programme's frames apart — and, on a render's last event, its log record. `true` continues, `false` stops:
`solve` answers the new `MasteringSolveStatus::Cancelled`, the range measurements `false`, and the same call made again
gives the same result. The overloads without a callback are unchanged, and the results are the same bits either way.

**C ABI v5.** `fc_master_set_progress (h, fc_progress_fn fn, void* context)` sets a handle's callback; it receives an
`fc_progress` (stage, pass, maxPasses, fraction, hasRecord, record) and returns 0 to stop, and the call then answers
`FC_ERR_CANCELLED` (15). On the wasm tier a handle without one calls `Module.onProgress(msg)` when that is a function;
`false` or an exception stops. No struct with a header grew. Checked on a real 5:21 stereo programme in node by
`tools/wasm/progress-check.mjs`.

**The render, too.** `DeliveredMastering::render` gained the same `progress` overload as its siblings, reporting
`Convert` then a new `Render` stage — no log record, since there is no search. `fc_master_render_delivered` takes the
handle's callback the same way `fc_master_solve_delivered` and `fc_master_measure_lra` already did (a new code,
`FC_PROGRESS_RENDER`; no struct grew) and cancels the same way. Checked in `tools/tests/MasterAbiProgressTests.cpp`.

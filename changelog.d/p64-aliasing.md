### `mastering` · `tools` — the loudness search refuses any input plane under an output plane; `fc_solution_log` refuses a `written` inside its records

- **Cross-channel aliasing is refused.** `TargetLoudnessSolver::solve` refused only `in[c] == out[c]`, so
  `out[0] = in[1]` was accepted: the render wrote channel 0's master where the next pass reads channel 1, and the
  call returned an ordinary verdict at a plausible gain over a master that is not the programme's — in the suite's
  witness the aliased call answered `TargetUnreachable` where the honest solve is `Solved`, at 12.2532 dB against
  12.3175, reporting −10.072 LUFS against −12.093, with a delivered master different in every one of 144 000
  frames. A search that happens to end after ONE render over such buffers was correct, and is refused too: whether
  it is correct would depend on how many passes it took, the reason `in == out` was already refused. The rule is
  now the one the C ABI and
  `DeliveredMastering` already applied: every plane non-null, and no input plane's bytes touching any output
  plane's, every pair, half-open. It is ONE predicate, `TargetLoudnessSolver::planesUsable`, which
  `DeliveredMastering` now forwards to instead of keeping its own copy. Still legal: one buffer feeding two input
  channels, planes edge to edge in one allocation, a call shorter than its buffers, buffers reused across calls —
  each pinned bit for bit against a solve on disjoint buffers. The direct C++ call now refuses what the facade
  refuses on the same memory (`MasterAbiTests`). **Also:** a single null plane beside good ones used to reach the
  renderer and crash; it is `InvalidRequest` now, from the same predicate.
- **`fc_solution_log`: `written` may not point into the `cap` records** (FC_ERR_SPAN), and is cleared only once
  every check is behind the call — the order `fc_master_flush` takes. It used to be cleared on entry, so a
  `written` inside the buffer took the count over a copied record on FC_OK and a zero into the buffer on a
  refusal. **A behaviour change for a caller that read `written` after a refused call:** a refusal on `out` (null,
  alignment, span) now leaves it as it was instead of zeroing it; `cap == 0` is FC_OK with a zero, as before. The
  refused span is the whole capacity, not the records a given log fills. Not a version: this change does not move
  `FC_MASTER_ABI_VERSION` (VERSIONING rule 1 — no struct grew and no entry point was added).
- Comments that said what the code does not: `renderTapped` claimed the chain's drain produces no gain reduction
  (it carries the release, and an expanding or upward mode acts on its silence; the windows exclude it, which is
  why the numbers were right); `fc_master.cpp` claimed a field retyped or inserted mid-struct is a build error (only
  where the change moves an offset or a size the pins read: an `int32_t` dropped into padding, `int32_t` to
  `uint32_t`, or a `double` narrowed to `float` before another `double` all build, and layout-check compares
  offsets, not types). And `fc_master_abi.h` said a count out-parameter is cleared FIRST, which no entry point with
  one does any more.

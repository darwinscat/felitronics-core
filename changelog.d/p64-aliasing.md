### `mastering` · `tools` — one rule for the planes a whole-programme operation reads and writes; `fc_solution_log` refuses a `written` inside its records

- **One rule, `mastering::planesUsable` (new `Planes.h`).** Tables and planes non-null; no input plane's bytes
  touching any output plane's, every pair, half-open, EACH SIDE AT ITS OWN LENGTH (`inFrames` for the input,
  `outFrames` for the output — a conversion's two differ, and judging both at one of them misses an overlap that lies
  only in the longer span or refuses planes that never meet); and no two output planes touching. The loudness search
  asks it with one length for both sides, `DeliveryConverter::convert` and `DeliveredMastering` with a conversion's
  two. Still legal: one buffer feeding two input channels, planes edge to edge in one allocation, a call shorter than
  its buffers, buffers reused across calls — each pinned bit for bit against buffers of their own.
- **The loudness search refused only `in[c] == out[c]`.** `out[0] = in[1]` was accepted: the render wrote channel 0's
  master where the next pass reads channel 1, and the call returned an ordinary verdict at a plausible gain over a
  master that is not the programme's — in the suite's witness the aliased call answered `TargetUnreachable` where the
  honest solve is `Solved`, at 12.2532 dB against 12.3175, reporting −10.072 LUFS against −12.093, with a delivered
  master different in every one of 144 000 frames. `out[0] = out[1]` was accepted too: `Solved` at 12.3041 dB against
  12.3175, with only the second channel's render left in the buffer where two channels were asked for. Both are
  `InvalidRequest` now, before a render. A search that happens to end after ONE render over cross-aliased buffers was correct, and is refused too:
  whether it is correct would depend on how many passes it took, the reason `in == out` was already refused. The
  direct C++ call now refuses what the facade refuses on the same memory (`MasterAbiTests`). **Also:** a single null
  plane beside good ones used to reach the renderer and crash; it is `InvalidRequest` now, from the same rule.
- **`DeliveryConverter::convert` checked its planes for null only.** It writes `out` at the delivery stride while still
  reading `in` at the source one, so an output plane over an input plane overwrote programme not yet read: at
  44.1 → 48 kHz `out[0] = in[1]` returned true with channel 1 wrong in 25 990 of 52 245 frames of the suite's
  witness; at 48 → 44.1 kHz it came out right only because there the writes lag the reads. Refused now, in both
  directions, before a sample is written — so is an overlap that exists only in the longer of the two spans, and two
  output planes on one buffer. **A behaviour change for a direct caller:** an identity conversion (equal rates) with
  `in[c] == out[c]` copied the bits correctly in place and is refused too — whether an overlap is safe would otherwise
  depend on the ratio; `DeliveredMastering` and the C ABI already refused it, and nothing in the tree calls it so.
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
  `uint32_t`, or a `double` narrowed to `float` before another `double` all build outside v4's type-pinned fields,
  and layout-check compares offsets, not types). And `fc_master_abi.h` said a count out-parameter is cleared FIRST,
  which no entry point with one does any more — including the sentence under `fc_solution_gr_trace` that contrasted
  it with `fc_solution_log`.

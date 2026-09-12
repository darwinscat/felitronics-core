### `tools` · `mastering` — the C ABI's compatibility rule (v2, v3): delivery at another rate, and parallel compression

`tools/fc_master_abi.h` states the rule every later ABI change follows, so that a version is a row in a table
rather than an event, and makes the first two bumps by it.

- **The rule** (VERSIONING, rules 1-8). One version for the whole ABI, moved by a struct that grows or an entry
  point that is added. A bump appends fields to the END of structs that begin with a header, each IN field with a
  default under which the call is the previous version bit for bit; nested structs and header-less array elements
  never grow. Tail padding is a named field. The size of every (struct, version) is ONE row of
  `FC_MASTER_STRUCT_SIZES`, published by `fc_master_sizeof(id, version)`. An IN struct of any known version is read
  over this build's defaults at the CALLER's size, and every bound and alias check uses that size; a version newer
  than the build is refused by decision. An OUT struct keeps the caller's header and gets exactly its size.
- **The frozen writers.** `fc_*_default(out)` cannot know the caller's size, so they are frozen at v1 (stamp v1,
  write v1's bytes); `fc_*_defaults(out)` write the caller's stamped version. A field set after a frozen writer
  lies past the stamp and is not read — the JS accessor refuses to write or read such a struct at all, and the
  CLI, the parity harness and the suites moved to the versioned writers. `fc_master_sizeof_params/config` are
  frozen at v1 too, so a v1 page fails on the version, not on a size.
- **v2 — `fc_master_config::deliveryRate`.** Non-zero makes a DELIVERING handle: SRC first, chain and solver at the
  delivery rate, equal rates included. New entry points `fc_master_delivered_frames`, `fc_master_render_delivered`
  and `fc_master_solve_delivered`; `fc_master_measure_lra` converts on such a handle; `process`/`flush`/`solve`
  refuse. The composition is one core class, `mastering::DeliveredMastering`, which the facade and the selftest's
  direct path both call. Budgets: `FC_NEED_SOLVE` / `FC_NEED_MEASURE_LRA` on a delivering handle include the
  converted programme and judge a range on the DELIVERED length; `render_delivered` allocates nothing.
- **v3 — `compressorMix`** (P60's field) at the end of `fc_master_params` and `fc_master_resolved`. Clamped by the
  core and read back; a non-finite mix is refused.
- **A non-finite input sample costs a delivered render what it costs a plain one.** `DeliveryConverter` now gates
  every input sample with the chain's own rule (NaN/inf -> 0, clamp +-1e6) before the conversion and counts it:
  a windowed sinc used to spread one NaN over its kernel, and the chain's gate then zeroed every delivered sample
  it reached. Converting a NaN is bit-identical to converting a zero in its place.
- **Refusals move nothing.** `DeliveredMastering` reaches every verdict — renderer, width, overlapping planes, and
  everything the solver would refuse before a pass (`TargetLoudnessSolver::admits`, factored out of `solve()` with
  no change of verdict) — before the converter writes a sample or a programme buffer is allocated.
- **Proof.** A render, a solve and a range through a v1 build (`179e1c5`) and this one are byte-identical on the
  same fixture; v1- and v2-stamped structs render the v3 caller's bits at the previous values; the delivered render
  and search are bit-identical to the core on down, up and equal rates; `need` equals the counted allocation;
  `tools/wasm/layout-check.mjs` holds every JS field offset against the compiler (wasm32 on that tier) and
  `layout-gate-test.mjs` the JS gates; `master-parity.mjs` compares the delivered paths wasm vs native at the
  delivered length. Product check, 30 s of music: 96 -> 44.1 kHz -14.003 LUFS / -1.044 dBTP, 192 -> 48 kHz
  -14.003 / -1.046, measured by `fcore_measure` at the delivery rate.
- **Site transition.** The rule cannot reach a page already shipped against v1: its loader requires `version === 1`.
  The site moves to v3 as one release of the worker, its layout copy and the module.


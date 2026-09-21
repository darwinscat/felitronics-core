<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools — every enum code and every field domain of the mastering ABI, published once and held against real calls

**No ABI change.** `FC_MASTER_ABI_VERSION` stays 8 and no struct moves: what follows is metadata a page reads
out of `tools/wasm/fc-master-layout.mjs`, beside the layouts that already live there.

**Ten enum lists, all of them.** `FC_FILTER_TYPE`, `FC_DETECTOR`, `FC_LINK_MODE`, `FC_COMP_MODE`, `FC_SHAPE`,
`FC_NOISE_SHAPING`, `FC_GR_STATISTIC` (with `Percentile`), `FC_GR_STAGE` and `FC_PROGRESS_STAGE` join
`FC_EQ_AXIS`, each an array whose INDEX is the code, and `FC_ENUMS` says which header enum each one mirrors.
`layout-check.mjs` reads the enumerators out of `tools/fc_master_abi.h` and compares every list entry by entry
— declaration order, value and letters — so a code added, renumbered or permuted in the header fails there
instead of silently renaming a menu. Consumers held five of these by hand and had just drifted on the sixth.

**`FC_DOMAINS` — 89 rows, one per input field** of `fc_master_config`, `fc_master_params` (the EQ bands, the
mono bass, the compressor, the clipper, the limiter and the dither) and `fc_loudness_request`. Each row carries
the unit, the admitted interval, what a value outside it does — `refuse` with the FC_ERR_* that names it,
`clamp`, `verdict` (`fc_master_solve` answers FC_OK and the SUMMARY carries `FC_SOLVE_INVALID_REQUEST`), `free`
where the code bounds nothing and only finiteness is checked, or `any` — what a non-finite value does, where
the applied value can be read back, and what the domain depends on. Bounds that move with the sample rate are
written as `0.49*sr` / `8000/sr` and evaluated by `domainBound (bound, sampleRate)`; the four that depend on
another FIELD say so rather than averaging (the mono-bass stage narrows the width to exactly 2; the oversample
factor stops at 16 with the limiter, 64 with the clipper alone and nowhere with neither; the compressor's
250 ms lookahead ceiling is that stage's own). **Silent clamps are marked**, because the difference decides an
interface: a refusal arrives on the call that made it, a clamp arrives as nothing at all. No field of this ABI
has a list-valued domain — `oversampleFactor` admits every integer in its interval and `dither.bits` refuses
nothing — and the table says so instead of inventing one. `layout-check.mjs` holds the table against the struct
layouts in both directions: a row must name a value field that exists, and EVERY input field of the three
structs must have a row, so a field added in a later version cannot arrive without a domain.

**`felitronics_master_domains_tests`** holds the table against the running ABI with a real C-ABI call per
bound, at two sample rates for every rate-dependent row. The table is the test's ARGUMENT and every probe is
derived from its own numbers, so a bound moved there without the code moving with it lands on the wrong side of
the real boundary. A clamp is pinned by READING BACK what was applied — `fc_master_resolved` where it publishes
one, `fc_master_eq_curve` for the EQ lane fields, the rendered audio for the rest — and asserting that the value
at the bound and beyond it are the same number while one step inside it is a different one; acceptance alone
would pass against a clamp anywhere at all. A refusal is probed twice — one step out and WELL out — because a
bound moved INWARDS still refuses one step past itself, and that near probe alone could not tell a delivery
floor of 8000 Hz from one of 22050. The six clamps nothing can pin are printed by name on every run.
1235 checks; a planted-mutation round over the table killed 56 of 59, and all three survivors are a row DEMOTED
to a weaker claim, which the suite's header names as what it does not catch. Eleven more planted violations —
an enum permuted in the header, a list shortened, a row naming a field that no longer exists — are all caught
by `layout-check.mjs`.

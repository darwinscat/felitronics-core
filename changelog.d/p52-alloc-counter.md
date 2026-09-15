<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### test_support · tools — one allocation counter, and it proves itself before the suite runs

"`process()` does not allocate" is the repo's oldest RT claim and it was being measured by sixty-one
private instruments, fifty of which could not see the allocations that matter most.

A suite proves that claim by replacing the global allocation functions and reading a delta. The set to
replace is **eight**, not two: C++17 routes any object whose `alignof` exceeds
`__STDCPP_DEFAULT_NEW_ALIGNMENT__` (16 on every desktop row, **8** on wasm32) through
`operator new(std::size_t, std::align_val_t)`, a different function. On `v0.33.0` the idiom had been copied
into 61 translation units; **11 installed the over-aligned form and 50 did not**, and both kinds said "no
heap allocation in the audio path" in the same words. `core::SeamAllocator<64>` — and through it every
convolver buffer — goes that way, and so does `eq::EqEngine`, whose `alignof` is 64 and whose own header
already said a two-form counter would not see its 331 KiB.

The hole was **latent**, and that is part of the finding rather than a softening of it: re-measuring with a
fixed instrument found nothing that had fallen through (see the last bullet). `modules/eq/tests/` is the
sharpest case — the suite the blindness would have hurt most makes zero over-aligned allocations of its
own, because it puts the engine on the stack.

- **`test_support/alloc_counter.h`** is now the only counter in the tree: all eight `new`s and all twelve
  `delete`s, the byte accounting (including MSVC's container padding and the plain-object exemption the
  `win` row found), a one-shot fault switch and a one-shot re-entry hook. The 61 copies are gone —
  **1390 lines deleted from the 61 suites against 348 added**, and one 363-line header in their place. Including the header IS installing it; there is no macro,
  because a macro can be forgotten and a forgotten one is silent, whereas including it twice in one
  executable is a duplicate-symbol link error on every linker in the matrix.
- **The instrument is proven in the run whose conclusions depend on it.** Before `main()`, the header asks
  the counter for one allocation through **each of the eight forms** and requires the number to move; any
  that goes unseen is named and the run aborts. Not a check the suite could fail to reach, and not in
  `report()` — two suites here total in a harness of their own and never call it. Eight probes rather than
  two because the distinction is not the obvious one: *deleting* a form is harmless (the standard's
  defaults forward `new[]`→`new`, nothrow→throwing, array-aligned→scalar-aligned, so the allocation is
  still counted — measured on libc++, libstdc++ and the UCRT alike), while a form that is **present and
  stops counting** is exactly the original defect in a new shape. One mutation per form: eight of eight
  red, each naming itself.
- **`tools/lint/check-alloc-counter.mjs`** keeps it that way: a private `operator new` anywhere under
  `modules/`, `tools/` or `test_support/` fails CI by name. Lexed rather than grepped — the words
  "operator new" are ordinary prose in this tree — with a 32-case self-test and two negative controls in
  the workflow, one for the lint's reach and one that deletes the two over-aligned lines from the shared
  counter and requires the binary to refuse to run. It anchors on the operator rather than on the return
  type and decides scope by brace depth, because a return-type matcher lets through every replacement whose
  signature is spelled differently — a newline after `void*`, a `[[nodiscard]]`, a trailing return type —
  and falsely flags a class's own allocator.
- **What the fixed instrument then found, stated as a number rather than a reassurance.** Every
  over-aligned allocation in the tree was traced: 1741 distinct call stacks across all 128 test binaries.
  They come from `prepare()`, `setIr()` and `Bank::build()` — `AlignedVector::assign` under a convolver's
  `prepare` accounts for 1305 of them — and **not one stack names `process`, `analyse` or `applyGain`**.
  A second stand gave each over-aligned allocation a weight of 10⁶, so one landing inside a measured
  region could not be absorbed by a tolerance; exactly two suites moved, and both were among the eleven
  whose counters already saw them. So the blindness was real and nothing had slipped through it.

**What the counter still cannot see is named in the header** rather than left to be discovered: it replaces
the C++ allocation functions, not `malloc`, and `pffft.c` calls one directly. Measured rather than assumed —
with `--wrap=malloc,calloc,realloc` on the gcc row, 200 `process()` calls through
`MatrixConvolver<PffftRealFft>` make **zero** C allocations (four in the whole run, all inside `prepare`).

No audio moved: the full verbose output of all 132 tests is identical to `v0.33.0` apart from wall-clock
timings, on arm64 macOS and on MSVC alike, and all 128 `N checks, M failures` lines match line for line.

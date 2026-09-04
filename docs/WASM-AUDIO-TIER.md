<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# The `wasm-audio` tier — gated, not aspirational

[`DSP-ARCHITECTURE.md`](DSP-ARCHITECTURE.md) §2 has always named a `wasm-audio` tier and always marked it
*aspirational*. [`P0-WASM-SPIKE.md`](P0-WASM-SPIKE.md) proved the core **can** live under Emscripten. This
document is what turned that into a gate: a `wasm-audio` CMake preset and a CI job that builds every default
module for wasm32 with exceptions and RTTI off and no pthreads, runs the whole self-test suite in node, and
audits every emitted artifact.

**Result: 74 of 74 suites pass in node, zero warnings, and all 78 emitted `.wasm` audit clean.** A second,
checked configuration (`SAFE_HEAP` + `ASSERTIONS=2` + stack-overflow checks) runs the same tree and is green
over its 73 non-perf suites — the perf suite is excluded there, and the artifact audit runs over the release
build. The one suite absent from both is `felitronics_core_rtstreams_tests`, which needs real threads and is
excluded by construction rather than left permanently red.

---

## What the gate proves — and what it does not

The task this came from asked for a gate where *"an attempt to add a thread, an allocation in `process()`, an
exception, or an incompatible dependency fails the build."* Three of those four were true as soon as the tier
built. The fourth — the thread — was **not**, and finding that out took a measurement rather than a reading
of the flags; it is true now, but only because of one extra link flag and with two limits worth naming.

| | caught how | when |
|---|---|---|
| an exception / RTTI | `-fno-exceptions -fno-rtti` | **build** (at instantiation — see below) |
| an incompatible dependency | compile / link error | **build** |
| a **thread** on a live code path | `--wrap=pthread_create` → `undefined symbol: __wrap_pthread_create` | **build** (link) |
| a thread in code never emitted | — | not caught |
| a `std::mutex` | a lint over module headers, not the artifact | build (textual) |
| an allocation on an **allocation-counted** path | the counter in `test_support/felitronics_test.h` | **run** |
| a denormal stall (law 8) | **not at all** | — |

**The thread claim was false until one link flag made it true, and that is worth spelling out.** Measured on
emsdk 6.0.9: `std::thread` with no `-pthread` **compiles and links cleanly** — Emscripten links pthread stubs
whose `pthread_create` returns `ENOTSUP` — and only aborts at runtime. `find_package(Threads REQUIRED)` also
*succeeds* there and yields an empty `Threads::Threads`, so a CMake dependency on threads does not force
`-pthread` either. And the artifact audit cannot see it: without `-pthread` the memory is unshared and the
glue is clean, because nothing threaded was linked — the audit proves *"not built with `-pthread`"*, which
is a different statement.

`LINKER:--wrap=pthread_create` closes it. libc++'s `std::thread` constructor is header-inline, so the
`pthread_create` reference lands in **our** object file, and `--wrap` renames it to an
`__wrap_pthread_create` nobody defines. Negative control, run against this tree: a `std::thread` added to
`core::Smoother::prepare` fails **8 links** with `wasm-ld: error: undefined symbol: __wrap_pthread_create`.
Zero false positives across all 78 links.

Two honest limits:

- **Dead code is not caught.** The first attempt at that negative control put the thread in an unused inline
  function and the build stayed green — nothing was emitted, so nothing referenced `pthread_create`. Moving
  it into a function tests actually call is what turned the build red.
- **`std::mutex` is not caught.** `mutex::lock` lives inside `libc++-noexcept.a` rather than our object, and
  emscripten's stub returns 0, so a lock is a silent no-op. The job therefore also greps module headers for
  `<thread>`, `<mutex>`, `<condition_variable>`, `<future>` and friends. That is a lint over text, not proof
  from the artifact, and it is labelled as such.

**The allocation gate is a genuine gain, but read the row carefully.** `felitronics_test.h` enforces its
allocation count only `#if defined(_LIBCPP_VERSION)` — on libstdc++ it degrades to informational, because a
global counter cannot separate our allocations from the standard library's there. Emscripten is libc++, so
on this tier the no-alloc rule is **enforced**, where the `ubuntu-latest` rows only record it. What it is
NOT is a blanket "no allocation in any `process()`": it fires only where a suite installs the counter and
asserts on it. `LaneDynamics::processBand`, for one, has no such assertion — an allocation added there stays
green. Instrumenting every RT entry point would make the broad claim true; until then the narrow one is what
this gate earns.

**Law 8 is untouched.** A denormal stall costs 10–100× on a CPU without FTZ; it is a runtime property of the
hardware and no build gate can see it. The `analysis::KWeightingFilter` finding from the P0 spike is still
open and this job does nothing about it.

---

## Policy and physics are separate on purpose

- **Policy** — `-fno-exceptions -fno-rtti -ffp-contract=off` — lives in the `wasm-audio` preset in
  `CMakePresets.json`. It is the tier's *contract*.
- **Physics** — stack size, heap growth, assertions — lives in `CMakeLists.txt` under `if(EMSCRIPTEN)`,
  because **any** Emscripten build of this tree needs it. Leaving it to a preset would mean a developer who
  typed `emcmake cmake` by hand got silently wrong answers.

All three physics values were measured, not guessed:

| flag | Emscripten default | why ours differs |
|---|---|---|
| `-sSTACK_SIZE=8388608` | 64 KB | the suite puts a ~200 KB `EqEngine` on the stack — the same fact the `/STACK:8388608` line for MSVC exists for. 8 MB matches the Unix default. |
| `-sALLOW_MEMORY_GROWTH=1` | fixed 16 MB heap | three suites legitimately need more and died `Aborted(OOM)`. |
| `-sASSERTIONS=1` | off in Release | without it a failure prints a bare `Aborted()` and nothing else. |

**A blown wasm stack does not reliably trap — it corrupts first.** That is the finding worth carrying
forward. At the 64 KB default, before any crash, `felitronics_analysis_tests` produced *wrong answers*:

```
FAIL: no frame before the hop elapses
FAIL: the hop reported is the delta that occurred (600), not the request (512)
RuntimeError: table index is out of bounds
```

A reported hop of 600 against a requested 512 is not a crash, it is a plausible-looking number. Two `eq`
suites then failed with `memory access out of bounds`. All five went green on stack size alone — no DSP was
touched. On a native target the same overflow would have hit a guard page.

---

## What building the tier found

Everything below was invisible on desktop and surfaced from one act: compiling the tree as `wasm-audio`.

1. **Seven more `throw`s than the record claimed.** The core's only `throw` was known
   (`core/Fft.h`, fixed separately). The tier also compiles the *tests*, and the counted aligned
   `operator new` in six convolution RT-safety suites plus the pffft null suite each throw `std::bad_alloc`.
   Same guard, same reason.

2. **The `-fno-exceptions` failure is INSTANTIATION-dependent** — the earlier write-up said the opposite.
   Measured: `#include <felitronics/core/Fft.h>` alone under `-fno-exceptions` is **clean**; naming
   `SeamAllocator<float>::allocate`, or any `AlignedVector`, is what errors. A gate that merely included
   headers would have proven nothing. This is why `felitronics_no_exceptions_probe` reuses the
   header-hygiene TU verbatim — that TU carries explicit template instantiations at its foot.

3. **A test that was UB on any 32-bit target.** `measurement/tests/AdversarialTests.cpp` asserted
   `nextPow2(size_t(1) << 63)`. Where `size_t` is 32 bits that shift is undefined, so the overflow guard the
   case exists to test went untested and the assertion failed on nonsense. It is now spelled from
   `std::numeric_limits<std::size_t>::digits`, which exercises the guard at **whatever** width it is
   compiled for — strictly more than the 64-bit-only version tested.

4. **Eight public headers were outside the strict-warning gate.** `HeaderHygiene.cpp` is supposed to include
   every public header; it reached 78 of 98, and 8 of the missing ones ship in the default build
   (`MultiResSpectrumPaneFast`, `NonUniformConvolver`, `LaneDynamics`, `BandBallistics`, `NoiseGate`,
   `RelativeLevel`, `MagnitudeCurve`, `XcorrAlign` — each landed with its feature branch and none was added
   here). Seven went in clean. The eighth did not, and it is the interesting one: `LaneDynamics.h:155`
   compared `rangeDb == 0.0`, which **gcc** rejects under `-Wfloat-equal` and **clang does not** (clang
   exempts a comparison against a literal). It is now `std::fabs (rangeDb) <= 0.0` — exactly equivalent on
   every input, and deliberately not the obvious `<= 0.0`, because `rangeDb`'s sign is the direction and
   that would have disengaged half the modes. The remaining 12 are the opt-in `fftpffft` / `nam` /
   `rigplayer` headers, gated with their options.

5. **The no-threads checker was blind in exactly the failing case.** `check-no-threads.mjs` parsed only the
   memory *section*. A `-pthread` build declares no memory section at all — the shared memory is
   **imported** — so it printed `memories declared: 0` and passed. It parses imports now, and both controls
   are recorded: a `-pthread` build fails on the wasm alone, a normal one passes.

6. **Its glue-string scan false-positived on every artifact.** A bare `pthread` substring matched
   `-sASSERTIONS`' own deprecation string (`"Module.pthreadMainPrefixURL option was removed"`), so all 78
   modules "failed". Replaced with tokens that actually discriminate — measured counts, threaded build vs a
   tier test vs the P0 web artifact: `PThread` 31/0/0, `ENVIRONMENT_IS_PTHREAD` 20/0/0,
   `_emscripten_thread` 22/0/0, `Atomics.` 5/0/0, `pthread_create` 3/0/0, `SharedArrayBuffer` 1/0/0,
   `new Worker` 1/0/0.

7. **A 32-bit target produced no new warnings at all.** Predicted to "slam" the build (`size_t` 32-bit,
   `long double` 128-bit quad); measured zero, under the full strict set, on both emsdk 6.0.9 and gcc 14.2.

8. **Three public class templates that no gate had ever entered.** Because a compiler skips an
   uninstantiated template member entirely, the strict-warning and `-fno-exceptions` gates only ever saw
   what `HeaderHygiene.cpp` explicitly instantiates. `convolution::NonUniformConvolver`,
   `convolution::MatrixConvolverNupc` and `analysis::RollingSpectrumTapT` were not on that list — including
   their headers does **nothing**, so a `throw` inside any of their members would have been invisible on
   every CI row. They are instantiated now.

---

## The CI job

One job, `wasm-audio-tier`, on **`ubuntu-latest`** — and the runner choice is load-bearing for its last step.
Apple's libm and musl (what Emscripten compiles in) return different doubles from `tan()` at 88.2 and
192 kHz, so an arm64-macOS reference would diverge from wasm for a reason that is not a regression. glibc
agrees with musl at every rate measured. See [`P0-WASM-SPIKE.md`](P0-WASM-SPIKE.md).

emsdk is installed from upstream rather than through a third-party action, so the only thing the job trusts
is emscripten itself, and **6.0.9** is the toolchain every number in the spike report was measured on. The
clone is pinned with `--branch 6.0.9`: a plain `--depth 1` clone takes emsdk *master* — the installer and the
release manifest — and pins only the SDK version argument, which is not the same thing.

The steps, and what each is for:

1. **build + `ctest`** — the tier itself.
2. **`check-no-threads.mjs --all`** — audits **every** emitted `.wasm`, not one representative artifact.
3. **the threading-header lint** — the `std::mutex` gap above.
4. **a second, CHECKED configuration** — `-sSAFE_HEAP=1 -sASSERTIONS=2 -sSTACK_OVERFLOW_CHECK=2`, ~15 s to
   build and ~17 s to run, all green but the perf suite (excluded: its assertions are *relative* wall-clock
   and SAFE_HEAP instruments every load and store, so it would measure the instrumentation). This is the
   configuration that names a stack overflow instead of letting it corrupt quietly.
5. **`tools/wasm/build.sh`** — CI had never compiled the P0 spike at all, which is why it broke unnoticed
   during the session that wrote it.
6. **native↔wasm NULL test** on a **generated** fixture — no private audio enters a public repo. The
   comparison surface is the pre-gate 400 ms block-energy vector plus the true-peak linear maximum, not the
   scalar LUFS: a gated scalar is discontinuous in its own inputs and cannot carry a bit-exactness claim.
   The **checked** artifact (`SAFE_HEAP` + `ASSERTIONS=2` + stack checks) is diffed as well — building it
   without running it would prove nothing, and the libsoxr precedent is precisely a clean build that died
   at runtime.

The fixture (`tools/wasm/make-fixture.mjs`) calls no transcendental — an integer xorshift PRNG and float32
arithmetic only — so its bytes are identical on every JS engine and node version (on a little-endian host: a
TypedArray uses the platform's byte order, and the native tool reads f32**le**, so the harness already
assumed that). It is 10 s of stereo at 48 kHz in
five deliberate sections: ordinary noise, **digital silence** (the absolute gate), louder noise, an fs/4
pattern sampled 45° off the crests, and a decaying burst that stops on a transient at the very last sample.

Two of those sections were designed against a measurement, not a hunch.

**The fs/4 section.** An obvious "alternating ±A" sits at *exactly* Nyquist, where the reconstruction maximum
**is** A — measured `tp == sp`, so the compared true peak came from the sample-peak floor and the polyphase
FIR was never under test at all. Sampling fs/4 45° off the crests puts every sample at A while the
reconstruction peaks near A·√2, which is what puts the filter's own output into the compared number.

**The ending.** The first version decayed a burst across the last second — and `0.98 × 0.9995^48000` is
`3.7e-11`, about −209 dBFS, so the "abrupt ending" was silence and tested nothing. It is now a full-scale hit
occupying the **last 16 samples**, louder than anything before it, so the file's true peak is decided there.
That is what makes it a test: disabling `finish()` in `fcore::Probe` drops the reported true peak from
**1.3936 to 1.0125** — 2.8 dB under-reported — and changes the compared line. With the drain in place, the
whole ending is exact.

Current result on that fixture: **97 block energies plus true peak plus sample peak, every bit equal**,
native arm64 (Apple clang) vs wasm32 — and identical again from the checked `SAFE_HEAP` build.

---

## Reproducing

```sh
source ~/emsdk/emsdk_env.sh
cmake --preset wasm-audio && cmake --build --preset wasm-audio && ctest --preset wasm-audio
node tools/wasm/check-no-threads.mjs --all build-wasm-audio

# the NULL test
./tools/wasm/build.sh
cmake --preset desktop && cmake --build --preset desktop --target fcore_measure
node tools/wasm/make-fixture.mjs fixture.f32 48000 2 10
./build/tools/fcore_measure blocks 48000 2 fixture.f32                        > native.txt
node tools/wasm/parity.mjs tools/wasm/build/fcprobe.node.js 48000 2 fixture.f32 > wasm.txt
diff native.txt wasm.txt          # empty output IS the acceptance criterion
```

`CMakePresets.json` also carries a `desktop` preset, because `-DCMAKE_BUILD_TYPE=Release` is not optional
here and is easy to leave off: without it `felitronics_spectrum_pane_perf_tests` fails, 19 s against 1 s — a
footgun already hit and written down once. It names `Release` in its build and test presets too, which is
what multi-config generators (Visual Studio) actually read; `CMAKE_BUILD_TYPE` alone means nothing there.
Its `binaryDir` is the conventional `build/`, so a build directory created before v0.25.0 will still carry
that era's cached `FELITRONICS_NAM_CORE_TAG`/`..._NAMZ_TAG` — `set(... CACHE ...)` never overwrites an
existing entry — and fail to compile against the old namz. `rm -rf build` once; it is not a defect of the
preset, the documented raw command hits it identically.

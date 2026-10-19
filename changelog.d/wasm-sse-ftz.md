### core — `FlushToZero.h` compiles in a wasm build with `-msse4.1`

A wasm consumer that wants Eigen's SSE path (NAM) builds with `-msimd128 -msse4.1`; emscripten then defines
`__SSE2__` and lowers SSE onto SIMD128. `FlushToZero.h` keyed its x86 path on `__SSE2__` and called `_mm_setcsr`,
which emscripten does not provide because WebAssembly has no MXCSR, so no file that included it compiled under those
flags. The header now tests for wasm first and takes the path a target without an FP-control register always took:
`ScopedFlushToZero` is a no-op, and law 8's software flush keeps denormals out, as it already did in every wasm
build. Native x86-64, arm64 and the existing wasm build preprocess to the same source as before, byte for byte; MSVC
defines neither wasm macro, so its branch is the old condition verbatim. The wasm-audio tier now builds every public
header under `-msimd128 -msse4.1` (`felitronics_wasm_sse_probe`) and runs `felitronics_core_wasm_sse_tests`: under
those flags subnormals survive inside `ScopedFlushToZero`, and `flushDenormal` / `flushPoison` still zero them.

### nam · rigplayer · poweramp — moved to felitronics-guitar-core (BREAKING)

`felitronics::nam`, `felitronics::rigplayer` and `felitronics::poweramp` now live in
[felitronics-guitar-core](https://github.com/darwinscat/felitronics-guitar-core), with their suites, their
third-party notices (NeuralAmpModelerCore, Eigen, nlohmann/json, namz) and their det-math manifest lines.
Target names, namespaces and headers are unchanged, so a consumer changes its CMake and nothing else:
make felitronics-core available first, then felitronics-guitar-core, and drop `FELITRONICS_WITH_NAM`.
Core now REFUSES `-DFELITRONICS_WITH_NAM=ON` with a message that says where the modules went, instead of
failing later on a missing target; a build directory that cached the option needs `-UFELITRONICS_WITH_NAM`.

What core keeps for that repository, and exports so it states the same rules rather than a copy of them:

- `cmake/FelitronicsPolicy.cmake` — FP contraction (law 10), the MSVC stack + UTF-8, the emscripten
  physics, the sanitizer option and the hygiene warning set, moved verbatim out of the top-level
  `CMakeLists.txt`. Core includes it itself; felitronics-guitar-core includes the same file.
- `felitronics::test_support` — the test harness as an INTERFACE target: `felitronics_test.h`, the one
  allocation counter, the true-peak oracles, and two new headers that carry the law-11 harnesses moved
  verbatim out of the mastering suites: `law11_call_contract.h` (P1–P7, `ADAPT`, `allProperties`) and
  `law11c_pause.h` (the pause-is-silence fixtures). `PowerAmpStage`'s address of both laws is tested
  there through these same headers.

CI loses the `nam` matrix axis (eight build-test rows become six) and the NAM builds in the ASan and gcc-14 rows.
`neural` (the inference seam) stays in core.

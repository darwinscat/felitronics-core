cmake_minimum_required(VERSION 3.22)
project(felitronics_core VERSION 0.25.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # so clangd / IDEs parse the headers at the right C++ standard

if(MSVC)
    add_compile_options(/utf-8)         # sources use UTF-8 in comments (—, φ, ≈, ·) — MSVC needs this flag
    add_link_options(/STACK:8388608)    # 8 MB stack (the Unix default). Tests put a ~200 KB EqEngine on the
                                        # stack; MSVC's 1 MB default overflows → SegFault. Match Unix.
endif()

# The `wasm-audio` tier's PHYSICS. Its POLICY (-fno-exceptions -fno-rtti) lives in the `wasm-audio` CMake
# preset; these three are needed by ANY emscripten build of this tree, so they are not left to a preset a
# developer might forget. All three were measured on emsdk 6.0.9, not guessed:
#   STACK_SIZE           emscripten defaults to 64 KB. The suite puts a ~200 KB EqEngine on the stack (the
#                        same fact the MSVC /STACK line above exists for), and in wasm a blown stack does NOT
#                        reliably trap: it first produced WRONG ANSWERS in the analyser tap tests (a reported
#                        hop of 600 against a requested 512) and only then "memory access out of bounds".
#                        8 MB = the Unix default the MSVC line already matches.
#   ALLOW_MEMORY_GROWTH  the default heap is a FIXED 16 MB; the NUPC stress suite and two measurement suites
#                        legitimately need more and died with Aborted(OOM).
#   ASSERTIONS           without it a failure prints a bare "Aborted()" and nothing else. A gate whose red
#                        light says nothing is a bad gate, and these binaries are never shipped. It is a
#                        cache variable because a directory-level link option lands AFTER
#                        CMAKE_EXE_LINKER_FLAGS, so a hard-coded =1 here would silently DOWNGRADE the
#                        checked configuration's -sASSERTIONS=2 (emscripten keeps the last spelling —
#                        verified in the generated link.txt). The checked CI step passes 2 through this.
set(FELITRONICS_WASM_ASSERTIONS "1" CACHE STRING
    "-sASSERTIONS level for emscripten binaries (2 in the checked configuration)")
if(EMSCRIPTEN)
    add_link_options(-sSTACK_SIZE=8388608 -sALLOW_MEMORY_GROWTH=1 -sASSERTIONS=${FELITRONICS_WASM_ASSERTIONS})
    # ...and the one flag that makes "adding a thread fails the build" TRUE rather than aspirational.
    # Without it emscripten links pthread STUBS whose pthread_create returns ENOTSUP, so std::thread
    # compiles, links, and only aborts at runtime — measured. --wrap turns the reference into an undefined
    # __wrap_pthread_create at link time: `wasm-ld: error: undefined symbol: __wrap_pthread_create`. It
    # works because libc++'s thread constructor is header-inline, so the reference lands in OUR object.
    # Measured zero false positives across all 78 links in this tree.
    # It does NOT catch std::mutex — those live inside libc++-noexcept.a and the stub silently succeeds —
    # so the CI job also lints module headers for <thread>/<mutex>/<condition_variable>/<future>.
    # Deliberately on for EVERY emscripten build of this tree, not only the preset: DSP-ARCHITECTURE.md 2
    # declares `wasm-audio` single-threaded, so a threaded wasm build of this core is not a thing we mean to
    # allow silently. If that ever changes it is one line, changed on purpose.
    add_link_options(LINKER:--wrap=pthread_create)
endif()

# ASan + UBSan for dev + a dedicated CI job. UBSan makes platform-divergent UB LOUD on any host: a
# convolution divide-by-zero on an unprepared engine was SILENT on Apple-Silicon (ARM SDIV /0 → 0) yet
# crashed on x86-64 (IDIV #DE) and would HardFault on bare-metal. `-fno-sanitize-recover=all` turns every
# UB into a hard failure (== a trap), so a future port to "who knows what" (Cortex-M …) stays honest.
# Flags are global so every module self-test is instrumented. Clang/GCC only (MSVC has no UBSan).
option(FELITRONICS_ENABLE_SANITIZERS "Build the self-tests with ASan + UBSan (hard-fail; Clang/GCC)" OFF)
if(FELITRONICS_ENABLE_SANITIZERS)
    if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU"))
        message(FATAL_ERROR "FELITRONICS_ENABLE_SANITIZERS requires Clang or GCC (MSVC has no UBSan).")
    endif()
    message(STATUS "felitronics-core: ASan + UBSan enabled (hard-fail on any UB)")
    add_compile_options(-g -fno-omit-frame-pointer
                        -fsanitize=address,undefined -fno-sanitize-recover=all)  # NOT -fsanitize=integer: the dither RNG wraps unsigned on purpose
    add_link_options(-fsanitize=address,undefined)
    add_compile_definitions(_GLIBCXX_ASSERTIONS)                                  # libstdc++ container bounds checks (inert on libc++; ASan covers OOB there)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-ftrivial-auto-var-init=pattern)                     # cheap uninitialised-stack tripwire
    endif()
endif()

# ----------------------------------------------------------------------------
# felitronics-core — shared, framework-agnostic, JUCE-free, real-time-safe C++ DSP core.
#
# Independent modules, each separately includable so a consumer pulls only what it needs (an EQ
# plugin must NOT drag in the neural runtime). Most modules are header-only INTERFACE targets;
# fftpffft and nam are compiled only when explicitly opted in. The full per-module map with key
# types is docs/CORE-OVERVIEW.md. See docs/DSP-ARCHITECTURE.md for the design rationale.
#
# Module          Namespace               Target                Deps
#   core          felitronics::core       felitronics::core     —            (Sample, sizes, Math, Smoother, flush)
#   analysis      felitronics::analysis   felitronics::analysis core         (SpectrumTap; FFT via the seam)
#   dynamics      felitronics::dynamics   felitronics::dynamics  core         (EnvelopeFollower · GainComputer · Compressor · NoiseGate · AutoLeveler · TransientShaper · ChannelLinker)
#   eq            felitronics::eq         felitronics::eq        core+analysis (matched biquads + Cytomic SVF + EqEngine)
#
# Consume (pinned by tag, exactly like JUCE) — JUCE-free, so this fetch never pulls JUCE:
#   include(FetchContent)
#   FetchContent_Declare(felitronics_core GIT_REPOSITORY <url> GIT_TAG <tag>)
#   FetchContent_MakeAvailable(felitronics_core)
#   target_link_libraries(app PRIVATE felitronics::eq felitronics::dynamics)
#
# A transitional `teq::core` alias + <teq/*.h> shims let products that still say teq:: repoint
# their FetchContent without touching source (see modules/eq).
# ----------------------------------------------------------------------------

# Build the JUCE-free module self-tests when this IS the top-level project, or a parent opts in.
# A consumer that FetchContent's felitronics-core does NOT inherit the tests.
option(FELITRONICS_BUILD_TESTS "Build felitronics-core self-tests" "${PROJECT_IS_TOP_LEVEL}")

# Optional COMPILED SIMD FFT backend (pffft). Default OFF, so a consumer
# that does not opt in stays header-only with zero external deps. When ON, felitronics::fftpffft is built and
# (with tests) its cross-backend NULL suite runs. See modules/fftpffft + THIRD_PARTY_NOTICES.md.
option(FELITRONICS_WITH_PFFFT "Build the optional compiled pffft SIMD FFT backend (felitronics::fftpffft)" OFF)

# Optional COMPILED NAM backend. Default OFF preserves the JUCE-free, external-dependency-free header-only
# path; ON fetches the pinned NeuralAmpModelerCore + namz sources and builds felitronics::nam.
option(FELITRONICS_WITH_NAM "Build the optional NeuralAmpModelerCore backend (felitronics::nam)" OFF)

if(FELITRONICS_BUILD_TESTS)
    enable_testing()
endif()

# Order matters: dependees before dependents (analysis -> core; dynamics -> core; eq -> core + analysis).
add_subdirectory(modules/core)
add_subdirectory(modules/analysis)
add_subdirectory(modules/dynamics)
add_subdirectory(modules/eq)
add_subdirectory(modules/fftpffft)   # OPTIONAL, compiled pffft SIMD backend — before convolution, whose CabConvolver may select it
add_subdirectory(modules/convolution)
add_subdirectory(modules/oversampling)
add_subdirectory(modules/limiter)
add_subdirectory(modules/neural)
add_subdirectory(modules/nam)   # OPTIONAL, compiled NAM backend — depends on neural + core; self-gates on FELITRONICS_WITH_NAM
add_subdirectory(modules/rigplayer)   # OPTIONAL, the .orbitrig pack player — header-only over nam + convolution + lineareq + eq; self-gates on FELITRONICS_WITH_NAM
add_subdirectory(modules/saturation)
add_subdirectory(modules/poweramp)   # tube power-amp trunk (SagEnvelope + TubeStage + PowerAmpStage, from OrbitCab) — needs eq + oversampling + saturation
add_subdirectory(modules/stereo)
add_subdirectory(modules/multiband)
add_subdirectory(modules/dynamiceq)
add_subdirectory(modules/deesser)
add_subdirectory(modules/dither)
add_subdirectory(modules/lineareq)
add_subdirectory(modules/measurement)   # OFFLINE IR-measurement math (sweep/Farina deconv/gate) — double, message-thread
add_subdirectory(modules/blend)   # OFFLINE mic-blend engine (mixer params + kernels + Σ(gain·strip) blend) — double coeffs, float buffers
add_subdirectory(modules/io)   # OFFLINE file I/O (WAV read/write, moved from OrbitCapture's oc/wav.hpp) — zero-dep, message-thread
add_subdirectory(modules/text)   # JUCE-free UTF-8 → ASCII romanization for filename slugs (German/ß + Cyrillic), from OrbitCapture NAM's Matrix.h — zero-dep, message-thread

# Dev tools (reference-validation CLIs) — top-level builds only.
if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(tools)
endif()

# ----------------------------------------------------------------------------
# Header hygiene — one TU that #includes EVERY public header (felitronics + the teq shims, key templates
# explicitly instantiated) under a strict, DOWNSTREAM-grade warning set with -Werror. Consumers compile
# our headers under juce::juce_recommended_warning_flags (-Wconversion, -Wfloat-equal, …); this target
# makes that a core gate so the class of "clean here, warns in TabbyEQ" bugs dies in CI, not downstream.
# Clang/GCC only (MSVC maps differently; the clang/gcc jobs are the gate).
if(FELITRONICS_BUILD_TESTS AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU"))
    add_executable(felitronics_header_hygiene test_support/HeaderHygiene.cpp)
    target_link_libraries(felitronics_header_hygiene PRIVATE
        felitronics::core felitronics::analysis felitronics::dynamics felitronics::eq
        felitronics::convolution felitronics::oversampling felitronics::limiter felitronics::neural
        felitronics::saturation felitronics::poweramp felitronics::stereo felitronics::multiband felitronics::dynamiceq
        felitronics::deesser felitronics::dither felitronics::lineareq felitronics::measurement felitronics::blend
        felitronics::io felitronics::text)
    if(FELITRONICS_WITH_NAM)
        target_link_libraries(felitronics_header_hygiene PRIVATE felitronics::nam felitronics::rigplayer)
    endif()
    target_compile_features(felitronics_header_hygiene PRIVATE cxx_std_20)
    target_compile_options(felitronics_header_hygiene PRIVATE
        -Werror -Wall -Wextra -Wpedantic
        -Wconversion -Wsign-conversion -Wfloat-equal -Wshadow -Wcast-align
        -Wzero-as-null-pointer-constant -Wunreachable-code -Wswitch-enum -Wextra-semi
        $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wshorten-64-to-32;-Wimplicit-float-conversion;-Wconditional-uninitialized;-Wconstant-conversion;-Wbool-conversion>)
    # When selected, felitronics::convolution propagates the uniform FELITRONICS_WITH_PFFFT definition
    # and backend target here, pushing MatrixConvolver<PffftRealFft> through the strict gate below.
    add_test(NAME felitronics_header_hygiene COMMAND felitronics_header_hygiene)
endif()

# ----------------------------------------------------------------------------
# -fno-exceptions / -fno-rtti PROBE — compile-only. Reuses the hygiene TU verbatim (every public header
# PLUS the explicit instantiations at its foot), because the failure this guards is INSTANTIATION-dependent:
# `throw` inside a class-template member body is diagnosed when the member is instantiated, not when the
# header is included. Measured on emsdk 6.0.9 + Apple clang 21 — including core/Fft.h alone is clean; naming
# SeamAllocator<float>::allocate (or any AlignedVector) is what errors. So a probe that merely #included
# would prove nothing.
# Why a NATIVE target when the wasm-audio job already builds this way: the wasm job is one toolchain
# (emscripten/clang). This runs on every existing CI row, so gcc and Apple clang are covered too, and it
# stays green on a machine with no emsdk. Building it IS the assertion — there is nothing to run.
# NOT linked against felitronics::nam / ::rigplayer: the NAM backend legitimately uses try/catch (its
# upstream sources throw), so it is outside the exception-free tier by design, not by oversight.
# Clang/GCC only — MSVC spells this /EHs-c- and warns rather than errors, so it is not a gate there.
if(FELITRONICS_BUILD_TESTS AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU"))
    add_library(felitronics_no_exceptions_probe OBJECT test_support/HeaderHygiene.cpp)
    target_link_libraries(felitronics_no_exceptions_probe PRIVATE
        felitronics::core felitronics::analysis felitronics::dynamics felitronics::eq
        felitronics::convolution felitronics::oversampling felitronics::limiter felitronics::neural
        felitronics::saturation felitronics::poweramp felitronics::stereo felitronics::multiband felitronics::dynamiceq
        felitronics::deesser felitronics::dither felitronics::lineareq felitronics::measurement felitronics::blend
        felitronics::io felitronics::text)
    target_compile_features(felitronics_no_exceptions_probe PRIVATE cxx_std_20)
    target_compile_options(felitronics_no_exceptions_probe PRIVATE -fno-exceptions -fno-rtti)
endif()

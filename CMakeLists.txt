cmake_minimum_required(VERSION 3.22)
project(felitronics_core VERSION 0.7.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # so clangd / IDEs parse the headers at the right C++ standard

if(MSVC)
    add_compile_options(/utf-8)         # sources use UTF-8 in comments (—, φ, ≈, ·) — MSVC needs this flag
    add_link_options(/STACK:8388608)    # 8 MB stack (the Unix default). Tests put a ~200 KB EqEngine on the
                                        # stack; MSVC's 1 MB default overflows → SegFault. Match Unix.
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
# plugin must NOT drag in the neural runtime). All 15 modules below are header-only INTERFACE targets
# (convolution / neural / lineareq included); the full per-module map with key types is
# docs/CORE-OVERVIEW.md. See docs/DSP-ARCHITECTURE.md for the design rationale.
#
# Module          Namespace               Target                Deps
#   core          felitronics::core       felitronics::core     —            (Sample, sizes, Math, Smoother, flush)
#   analysis      felitronics::analysis   felitronics::analysis core         (SpectrumTap; FFT via the seam)
#   dynamics      felitronics::dynamics   felitronics::dynamics  core         (EnvelopeFollower · GainComputer · Compressor · NoiseGate · TransientShaper · ChannelLinker)
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

# Optional COMPILED SIMD FFT backend (pffft) — the one non-header-only module. Default OFF, so a consumer
# that does not opt in stays header-only with zero external deps. When ON, felitronics::fftpffft is built and
# (with tests) its cross-backend NULL suite runs. See modules/fftpffft + THIRD_PARTY_NOTICES.md.
option(FELITRONICS_WITH_PFFFT "Build the optional compiled pffft SIMD FFT backend (felitronics::fftpffft)" OFF)

if(FELITRONICS_BUILD_TESTS)
    enable_testing()
endif()

# Order matters: dependees before dependents (analysis -> core; dynamics -> core; eq -> core + analysis).
add_subdirectory(modules/core)
add_subdirectory(modules/analysis)
add_subdirectory(modules/dynamics)
add_subdirectory(modules/eq)
add_subdirectory(modules/convolution)
add_subdirectory(modules/oversampling)
add_subdirectory(modules/limiter)
add_subdirectory(modules/neural)
add_subdirectory(modules/saturation)
add_subdirectory(modules/stereo)
add_subdirectory(modules/multiband)
add_subdirectory(modules/dynamiceq)
add_subdirectory(modules/deesser)
add_subdirectory(modules/dither)
add_subdirectory(modules/lineareq)
add_subdirectory(modules/measurement)   # OFFLINE IR-measurement math (sweep/Farina deconv/gate) — double, message-thread
add_subdirectory(modules/blend)   # OFFLINE mic-blend engine (mixer params + kernels + Σ(gain·strip) blend) — double coeffs, float buffers
add_subdirectory(modules/fftpffft)   # OPTIONAL, compiled (pffft SIMD backend) — self-gates on FELITRONICS_WITH_PFFFT

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
        felitronics::saturation felitronics::stereo felitronics::multiband felitronics::dynamiceq
        felitronics::deesser felitronics::dither felitronics::lineareq felitronics::measurement felitronics::blend)
    target_compile_features(felitronics_header_hygiene PRIVATE cxx_std_20)
    target_compile_options(felitronics_header_hygiene PRIVATE
        -Werror -Wall -Wextra -Wpedantic
        -Wconversion -Wsign-conversion -Wfloat-equal -Wshadow -Wcast-align
        -Wzero-as-null-pointer-constant -Wunreachable-code -Wswitch-enum -Wextra-semi
        $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wshorten-64-to-32;-Wimplicit-float-conversion;-Wconditional-uninitialized;-Wconstant-conversion;-Wbool-conversion>)
    if(FELITRONICS_WITH_PFFFT)   # push MatrixConvolver<PffftRealFft> through the strict gate (PR2 C1 asserts + row stride)
        target_compile_definitions(felitronics_header_hygiene PRIVATE FELITRONICS_WITH_PFFFT)
        target_link_libraries(felitronics_header_hygiene PRIVATE felitronics::fftpffft)
    endif()
    add_test(NAME felitronics_header_hygiene COMMAND felitronics_header_hygiene)
endif()

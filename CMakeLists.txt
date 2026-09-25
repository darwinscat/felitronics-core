cmake_minimum_required(VERSION 3.22)
project(felitronics_core VERSION 0.51.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # so clangd / IDEs parse the headers at the right C++ standard

# The build policy (FP contraction, MSVC stack + UTF-8, emscripten physics, sanitizers, the hygiene
# warning set) lives in one file so a satellite repository can state the SAME flags — see its header.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FelitronicsPolicy.cmake)

# ----------------------------------------------------------------------------
# felitronics-core — shared, framework-agnostic, JUCE-free, real-time-safe C++ DSP core.
#
# Independent modules, each separately includable so a consumer pulls only what it needs (an EQ
# plugin must NOT drag in the neural runtime). Most modules are header-only INTERFACE targets;
# fftpffft is compiled only when explicitly opted in. The full per-module map with key
# types is docs/CORE-OVERVIEW.md. See docs/DSP-ARCHITECTURE.md for the design rationale.
#
# Module          Namespace               Target                Deps
#   core          felitronics::core       felitronics::core     —            (Sample, sizes, Math, Smoother, flush)
#   analysis      felitronics::analysis   felitronics::analysis core+oversampling (SpectrumTap; FFT via the seam; the reference true peak)
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

# The NAM backend, the .orbitrig pack player and the tube power-amp trunk (felitronics::nam, ::rigplayer,
# ::poweramp) moved to felitronics-guitar-core. A build that still asks for them is told so, instead of
# failing later on a target that no longer exists here.
if(FELITRONICS_WITH_NAM)
    message(FATAL_ERROR
        "FELITRONICS_WITH_NAM: felitronics::nam, ::rigplayer and ::poweramp moved to felitronics-guitar-core "
        "(https://github.com/darwinscat/felitronics-guitar-core). Make felitronics-core available first, then "
        "that repository, and drop FELITRONICS_WITH_NAM from your CMake and your cache (-UFELITRONICS_WITH_NAM).")
endif()

if(FELITRONICS_BUILD_TESTS)
    enable_testing()
endif()

# Order matters: dependees before dependents (oversampling -> core; analysis -> core + oversampling; dynamics -> core;
# eq -> core + analysis).
add_subdirectory(modules/core)
add_subdirectory(modules/oversampling)
add_subdirectory(modules/analysis)
add_subdirectory(modules/dynamics)
add_subdirectory(modules/eq)
add_subdirectory(modules/fftpffft)   # OPTIONAL, compiled pffft SIMD backend — before convolution, whose CabConvolver may select it
add_subdirectory(modules/convolution)
add_subdirectory(modules/limiter)
add_subdirectory(modules/neural)
add_subdirectory(modules/saturation)
add_subdirectory(modules/stereo)
add_subdirectory(modules/multiband)
add_subdirectory(modules/dynamiceq)
add_subdirectory(modules/deesser)
add_subdirectory(modules/dither)
add_subdirectory(modules/mastering)   # the mastering CHAIN — after every stage it composes (eq, dynamics, saturation, limiter, stereo, dither)
add_subdirectory(modules/lineareq)
add_subdirectory(modules/measurement)   # OFFLINE IR-measurement math (sweep/Farina deconv/gate) — double, message-thread
add_subdirectory(modules/blend)   # OFFLINE mic-blend engine (mixer params + kernels + Σ(gain·strip) blend) — double coeffs, float buffers
add_subdirectory(modules/io)   # OFFLINE file I/O (WAV read/write, moved from OrbitCapture's oc/wav.hpp) — zero-dep, message-thread
add_subdirectory(modules/text)   # JUCE-free UTF-8 → ASCII romanization for filename slugs (German/ß + Cyrillic), from OrbitCapture NAM's Matrix.h — zero-dep, message-thread

# felitronics::test_support — the house test harness (felitronics_test.h, the ONE allocation counter, the
# true-peak oracles, the law-11 property harnesses) as a target. Defined whether or not core's own tests
# are built, so a satellite repository's suites link it rather than reach into this tree by path — and
# get the same counter, not a second one.
add_library(felitronics_test_support INTERFACE)
add_library(felitronics::test_support ALIAS felitronics_test_support)
target_include_directories(felitronics_test_support INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/test_support)
target_link_libraries(felitronics_test_support INTERFACE felitronics::core)
target_compile_features(felitronics_test_support INTERFACE cxx_std_20)

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
        felitronics::deesser felitronics::dither felitronics::mastering felitronics::lineareq felitronics::measurement felitronics::blend
        felitronics::io felitronics::text)
    target_compile_features(felitronics_header_hygiene PRIVATE cxx_std_20)
    target_compile_options(felitronics_header_hygiene PRIVATE ${FELITRONICS_HYGIENE_WARNINGS})
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
# Clang/GCC only — MSVC spells this /EHs-c- and warns rather than errors, so it is not a gate there.
if(FELITRONICS_BUILD_TESTS AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU"))
    add_library(felitronics_no_exceptions_probe OBJECT test_support/HeaderHygiene.cpp)
    target_link_libraries(felitronics_no_exceptions_probe PRIVATE
        felitronics::core felitronics::analysis felitronics::dynamics felitronics::eq
        felitronics::convolution felitronics::oversampling felitronics::limiter felitronics::neural
        felitronics::saturation felitronics::stereo felitronics::multiband felitronics::dynamiceq
        felitronics::deesser felitronics::dither felitronics::mastering felitronics::lineareq felitronics::measurement felitronics::blend
        felitronics::io felitronics::text)
    target_compile_features(felitronics_no_exceptions_probe PRIVATE cxx_std_20)
    target_compile_options(felitronics_no_exceptions_probe PRIVATE -fno-exceptions -fno-rtti)
endif()

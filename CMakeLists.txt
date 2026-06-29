cmake_minimum_required(VERSION 3.22)
project(felitronics_core VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # so clangd / IDEs parse the headers at the right C++ standard

# ----------------------------------------------------------------------------
# felitronics-core — shared, framework-agnostic, JUCE-free, real-time-safe C++ DSP core.
#
# Independent modules, each separately includable so a consumer pulls only what it needs (an EQ
# plugin must NOT drag in the neural runtime). The LIGHT modules below are header-only INTERFACE
# targets; the HEAVY / platform-specific ones (convolution, neural, FFT backends) land later as
# compiled STATIC/OBJECT targets with narrow public headers. See docs/DSP-ARCHITECTURE.md.
#
# Module          Namespace               Target                Deps
#   core          felitronics::core       felitronics::core     —            (Sample, sizes, Math, Smoother, flush)
#   analysis      felitronics::analysis   felitronics::analysis core         (SpectrumTap; FFT via the seam)
#   dynamics      felitronics::dynamics   felitronics::dynamics  core         (EnvelopeFollower + GainComputer)
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

# Dev tools (reference-validation CLIs) — top-level builds only.
if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(tools)
endif()

# felitronics::eq — matched biquads (Vicanek) + Cytomic SVF + EqBand + EqEngine, migrated verbatim
# from TabbyEQ's teq/ into the felitronics::eq namespace. Header-only; depends on core
# (Math/Smoother/sizes/flush) + analysis (SpectrumTap). The kMaxChannels SSOT lives in felitronics::core.
add_library(felitronics_eq INTERFACE)
add_library(felitronics::eq ALIAS felitronics_eq)
target_include_directories(felitronics_eq INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_eq INTERFACE felitronics::core felitronics::analysis)
target_compile_features(felitronics_eq INTERFACE cxx_std_20)

# --- Transitional teq:: compat ----------------------------------------------------------------------
# A REAL target named `teq_core` (NOT just an alias): OrbitCab reads the *target* `teq_core` via
# get_target_property(... teq_core INTERFACE_INCLUDE_DIRECTORIES) to mark teq's headers SYSTEM, so an
# alias alone would break its configure. teq_core carries felitronics::eq + the <teq/*.h> re-export
# shims (same include dir) → a consumer repoints its FetchContent URL/tag and keeps
# `#include <teq/EqEngine.h>` / `teq::EqEngine` with zero source change.
add_library(teq_core INTERFACE)
add_library(teq::core ALIAS teq_core)
target_link_libraries(teq_core INTERFACE felitronics::eq)
# Expose every shim-transitive include dir on teq_core's OWN INTERFACE_INCLUDE_DIRECTORIES, so a consumer
# that marks teq's headers SYSTEM by reading that property (OrbitCab's get_target_property hack) also
# covers the felitronics/core + felitronics/analysis headers the <teq/*.h> shims pull in.
target_include_directories(teq_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${felitronics_core_SOURCE_DIR}/modules/core/include
    ${felitronics_core_SOURCE_DIR}/modules/analysis/include)

if(FELITRONICS_BUILD_TESTS)
    # The TabbyEQ teq test-suite (227 checks), migrated VERBATIM and run against the teq:: shims — a
    # pass proves the felitronics::eq re-home is behaviour-identical to the proven teq:: core.
    add_executable(felitronics_eq_tests
        tests/TestMain.cpp
        tests/MatchedBiquadTests.cpp
        tests/EqEngineTests.cpp)
    target_link_libraries(felitronics_eq_tests PRIVATE felitronics::eq)
    target_include_directories(felitronics_eq_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_compile_features(felitronics_eq_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_eq_tests COMMAND felitronics_eq_tests)
endif()

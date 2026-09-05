# felitronics::mastering — the mastering CHAIN, not new DSP: gain -> EQ -> [mono-bass] -> compressor
# -> [clipper] -> gain -> true-peak limiter -> dither, plus the offline render wrapper. It is the one
# module that deliberately depends on many others (see docs/DSP-ARCHITECTURE.md §4): the composite is
# the unit under test here, while the VOICING — preset tables, target curves, defaults — stays in the
# product. Header-only.
add_library(felitronics_mastering INTERFACE)
add_library(felitronics::mastering ALIAS felitronics_mastering)
target_include_directories(felitronics_mastering INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_mastering INTERFACE
    felitronics::core
    felitronics::eq
    felitronics::dynamics
    felitronics::saturation
    felitronics::limiter
    felitronics::stereo
    felitronics::dither)
target_compile_features(felitronics_mastering INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The chain's own acceptance: block-invariance on an adversarial partition, the bypass null, the
    # latency identity, and the sequences (channel-count refusal, reset, failed prepare, oversized
    # block, non-finite input).
    add_executable(felitronics_mastering_tests tests/MasteringChainTests.cpp)
    target_link_libraries(felitronics_mastering_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_mastering_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_mastering_tests COMMAND felitronics_mastering_tests)

    # Kept separate because it is a different KIND of test: it re-derives the chain's promises from
    # outside — a hand-built reference chain, an independently computed alignment oracle, and the
    # measured passband table — so a defect that is consistent with the chain's own arithmetic still
    # fails here.
    add_executable(felitronics_mastering_oracle_tests tests/MasteringOracleTests.cpp)
    target_link_libraries(felitronics_mastering_oracle_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_oracle_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_mastering_oracle_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_mastering_oracle_tests COMMAND felitronics_mastering_oracle_tests)
endif()

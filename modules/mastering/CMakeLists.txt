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
    felitronics::dither
    # The loudness solver measures what it renders, so the module now owns a METER dependency as well as
    # a chain of processors. It is the same argument the chain itself makes: the composite is the unit,
    # and a solver that could not measure would be a solver a product had to finish.
    felitronics::analysis)
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

    # LAW 11 — the caller's contract, as PROPERTIES over the whole accepted and the whole rejected
    # domain of every block-level entry point in the core. It lives here rather than in a module
    # because the law is cross-module by construction: the census that produced it found three answers
    # for the width and four for the length, and a single-module suite is exactly what let that happen.
    # Links every module that has such an entry point (mastering already pulls most of them in).
    add_executable(felitronics_call_contract_tests tests/CallContractTests.cpp)
    target_link_libraries(felitronics_call_contract_tests PRIVATE
        felitronics::mastering felitronics::multiband felitronics::deesser felitronics::dynamiceq
        felitronics::poweramp felitronics::lineareq felitronics::convolution felitronics::analysis)
    target_include_directories(felitronics_call_contract_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_call_contract_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_call_contract_tests COMMAND felitronics_call_contract_tests)

    # Law 11c ("a pause is silence") spans seven addresses in six modules plus one composite, and the law
    # is ONE claim — split per module, each stage would get to answer it its own way. Lives here for the
    # same reason CallContractTests does: this is the module that already links everything law 11 touches.
    add_executable(felitronics_pause_is_silence_tests tests/PauseIsSilenceTests.cpp)
    target_link_libraries(felitronics_pause_is_silence_tests PRIVATE
        felitronics::mastering felitronics::multiband felitronics::deesser felitronics::dynamiceq
        felitronics::poweramp felitronics::eq)
    target_include_directories(felitronics_pause_is_silence_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_pause_is_silence_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_pause_is_silence_tests COMMAND felitronics_pause_is_silence_tests)

    # P7 — the target-loudness solver, its constraints and its statistics. Separate from the chain's own
    # suite because it tests a different KIND of claim: not "the chain computes this", but "the search
    # lands here, in this many renders, and refuses by name when it cannot".
    add_executable(felitronics_loudness_solver_tests tests/LoudnessSolverTests.cpp)
    target_link_libraries(felitronics_loudness_solver_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_loudness_solver_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_loudness_solver_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_loudness_solver_tests COMMAND felitronics_loudness_solver_tests)
endif()

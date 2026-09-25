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
    # The EQ stage's dynamics producer: each band's per-lane delta seam is driven by
    # `dynamiceq::LaneDynamics`, which `eq` itself deliberately does not know about.
    felitronics::dynamiceq
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

    # The EQ stage's dynamics, as a NULL against the producer built outside this chain: the oracle is
    # `eq::EqBand` + `dynamiceq::LaneDynamics` on the same signal, and the unarmed case is
    # `eq::EqEngine::process` over the same quantum.
    add_executable(felitronics_mastering_dyneq_tests tests/DynamicEqInChainTests.cpp)
    target_link_libraries(felitronics_mastering_dyneq_tests PRIVATE felitronics::mastering felitronics::dynamiceq)
    target_include_directories(felitronics_mastering_dyneq_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_mastering_dyneq_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_mastering_dyneq_tests COMMAND felitronics_mastering_dyneq_tests)

    # Kept separate because it is a different KIND of test: it re-derives the chain's promises from
    # outside — a hand-built reference chain, an independently computed alignment oracle, and the
    # measured passband table — so a defect that is consistent with the chain's own arithmetic still
    # fails here.
    # P57 step 2 — SRC first: the whole-programme delivery conversion both render paths call.
    add_executable(felitronics_mastering_deliveryconverter_tests tests/DeliveryConverterTests.cpp)
    target_link_libraries(felitronics_mastering_deliveryconverter_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_deliveryconverter_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    add_test(NAME felitronics_mastering_deliveryconverter_tests COMMAND felitronics_mastering_deliveryconverter_tests)

    # P57b — the three delivered operations (render, solve, range) as ONE composition both render paths call.
    add_executable(felitronics_mastering_deliveredmastering_tests tests/DeliveredMasteringTests.cpp)
    target_link_libraries(felitronics_mastering_deliveredmastering_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_deliveredmastering_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    add_test(NAME felitronics_mastering_deliveredmastering_tests COMMAND felitronics_mastering_deliveredmastering_tests)

    add_executable(felitronics_mastering_oracle_tests tests/MasteringOracleTests.cpp)
    target_link_libraries(felitronics_mastering_oracle_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_oracle_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_mastering_oracle_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_mastering_oracle_tests COMMAND felitronics_mastering_oracle_tests)

    # P60 — `compressorMix`, parallel compression. Its own binary because its three gates are claims about
    # ONE field that the chain's suite would otherwise scatter: mix 1 is the compressor alone, mix 0 is the
    # input found in the input, and between them the blend, within a float ulp of exact and in the double
    # arithmetic the header states.
    add_executable(felitronics_compressor_mix_tests tests/CompressorMixTests.cpp)
    target_link_libraries(felitronics_compressor_mix_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_compressor_mix_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_compressor_mix_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_compressor_mix_tests COMMAND felitronics_compressor_mix_tests)

    # P62 — how far apart the core's two true-peak meters read, on delivered masters at every delivery rate. The
    # only owner of that number: headers that mention the gap point here instead of quoting it.
    add_executable(felitronics_truepeak_instrument_gap_tests tests/TruePeakInstrumentGapTests.cpp)
    target_link_libraries(felitronics_truepeak_instrument_gap_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_truepeak_instrument_gap_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_truepeak_instrument_gap_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_truepeak_instrument_gap_tests COMMAND felitronics_truepeak_instrument_gap_tests)

    # P62 — every pair of source and delivery rates delivers under its promise as the CERTIFICATE reads it, and the
    # solver's report is that certificate bit for bit.
    add_executable(felitronics_delivered_ceiling_tests tests/DeliveredCeilingTests.cpp)
    target_link_libraries(felitronics_delivered_ceiling_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_delivered_ceiling_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_delivered_ceiling_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_delivered_ceiling_tests COMMAND felitronics_delivered_ceiling_tests)

    # P7 — the target-loudness solver, its constraints and its statistics. Separate from the chain's own
    # suite because it tests a different KIND of claim: not "the chain computes this", but "the search
    # lands here, in this many renders, and refuses by name when it cannot".
    add_executable(felitronics_band_gr_tests tests/BandGrTests.cpp)
    target_link_libraries(felitronics_band_gr_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_band_gr_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_band_gr_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_band_gr_tests COMMAND felitronics_band_gr_tests)

    add_executable(felitronics_loudness_solver_tests tests/LoudnessSolverTests.cpp)
    target_link_libraries(felitronics_loudness_solver_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_loudness_solver_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_loudness_solver_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_loudness_solver_tests COMMAND felitronics_loudness_solver_tests)

    add_executable(felitronics_mastering_progress_tests tests/ProgressTests.cpp)
    target_link_libraries(felitronics_mastering_progress_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_progress_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_mastering_progress_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_mastering_progress_tests COMMAND felitronics_mastering_progress_tests)

    # Law 11b for OfflineRenderer and the chain — the named case the law-11 census carried for them, through
    # the same harness (test_support/law11_call_contract.h); the census itself is modules/laws.
    add_executable(felitronics_mastering_call_contract_tests tests/MasteringCallContractTests.cpp)
    target_link_libraries(felitronics_mastering_call_contract_tests PRIVATE felitronics::mastering)
    target_include_directories(felitronics_mastering_call_contract_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_mastering_call_contract_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_mastering_call_contract_tests COMMAND felitronics_mastering_call_contract_tests)
endif()

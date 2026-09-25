# modules/laws — NOT a module: no library, no headers, no namespace. It holds the two suites that test the
# core's LAWS across modules rather than any one module, and it is a directory under modules/ on purpose:
# tools/ci/scope.py derives the selection of a pull request from the tree (a module's users are whoever
# #include its headers, tests included), so these suites are selected whenever any module they drive
# changes — with no hand-kept list, and a directory outside modules/ would be selected only by a change to
# itself. Tests only; nothing here is built for a consumer.
if(FELITRONICS_BUILD_TESTS)
    # LAW 11 — the caller's contract, as PROPERTIES over the whole accepted and the whole rejected
    # domain of every block-level entry point in the core. It lives outside any one module because the law
    # is cross-module by construction: the census that produced it found three answers for the width and
    # four for the length, and a single-module suite is exactly what let that happen. Links every module
    # that has such an entry point.
    add_executable(felitronics_call_contract_tests tests/CallContractTests.cpp)
    target_link_libraries(felitronics_call_contract_tests PRIVATE
        felitronics::analysis felitronics::convolution felitronics::deesser felitronics::dither
        felitronics::dynamiceq felitronics::dynamics felitronics::eq felitronics::limiter
        felitronics::lineareq felitronics::multiband felitronics::saturation felitronics::stereo)
    target_include_directories(felitronics_call_contract_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_call_contract_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_call_contract_tests COMMAND felitronics_call_contract_tests)

    # Law 11c ("a pause is silence") spans six addresses in five modules plus one composite, and the law
    # is ONE claim — split per module, each stage would get to answer it its own way. Here for the same
    # reason CallContractTests is. The seventh address, poweramp::PowerAmpStage, moved with its module to
    # felitronics-guitar-core and is tested there through the same fixture helpers
    # (test_support/law11c_pause.h), not a copy of them.
    add_executable(felitronics_pause_is_silence_tests tests/PauseIsSilenceTests.cpp)
    target_link_libraries(felitronics_pause_is_silence_tests PRIVATE
        felitronics::deesser felitronics::dynamiceq felitronics::dynamics felitronics::eq
        felitronics::multiband)
    target_include_directories(felitronics_pause_is_silence_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_pause_is_silence_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_pause_is_silence_tests COMMAND felitronics_pause_is_silence_tests)
endif()

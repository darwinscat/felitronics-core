# felitronics::dehum — mains hum removal (50/60 Hz and harmonics) by a phase-coherent synchronous canceller,
# a constant-width notch comb. Header-only, JUCE-free, depends on felitronics::core only.
add_library(felitronics_dehum INTERFACE)
add_library(felitronics::dehum ALIAS felitronics_dehum)
target_include_directories(felitronics_dehum INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_dehum INTERFACE felitronics::core)
target_compile_features(felitronics_dehum INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The canceller as a filter: null depth, notch shape against the exact transfer, passband, ring.
    add_executable(felitronics_dehum_tests tests/HumRemoverTests.cpp)
    target_link_libraries(felitronics_dehum_tests PRIVATE felitronics::dehum)
    target_include_directories(felitronics_dehum_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_dehum_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_dehum_tests COMMAND felitronics_dehum_tests)

    # The canceller as a stage: laws 8/8a/11a/11c, digital black, poison, mode ramps, no allocation.
    add_executable(felitronics_dehum_contract_tests tests/HumRemoverContractTests.cpp)
    target_link_libraries(felitronics_dehum_contract_tests PRIVATE felitronics::dehum)
    target_include_directories(felitronics_dehum_contract_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_dehum_contract_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_dehum_contract_tests COMMAND felitronics_dehum_contract_tests)
endif()

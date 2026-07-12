# felitronics::saturation — oversampled soft-saturation (WaveShaper transfer curves + the Saturator that
# anti-alias-wraps them). Header-only, JUCE-free. Depends on core (Math) + oversampling (anti-alias).
add_library(felitronics_saturation INTERFACE)
add_library(felitronics::saturation ALIAS felitronics_saturation)
target_include_directories(felitronics_saturation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_saturation INTERFACE felitronics::core felitronics::oversampling)
target_compile_features(felitronics_saturation INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_saturation_tests tests/SaturationTests.cpp)
    target_link_libraries(felitronics_saturation_tests PRIVATE felitronics::saturation)
    target_include_directories(felitronics_saturation_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_saturation_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_saturation_tests COMMAND felitronics_saturation_tests)

    add_executable(felitronics_saturation_theory_tests tests/SaturatorTheoryTests.cpp)   # theory-first falsification: gate/symmetry/horizon/chunking contracts
    target_link_libraries(felitronics_saturation_theory_tests PRIVATE felitronics::saturation)
    target_include_directories(felitronics_saturation_theory_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_saturation_theory_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_saturation_theory_tests COMMAND felitronics_saturation_theory_tests)
endif()

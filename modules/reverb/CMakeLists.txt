# felitronics::reverb — Freeverb (the published algorithm, matched to the behaviour the amps were voiced on) and
# AmpReverb, the amp plugins' six-character tail. Header-only, JUCE-free; depends on core and eq (the wet HPF).
add_library(felitronics_reverb INTERFACE)
add_library(felitronics::reverb ALIAS felitronics_reverb)
target_include_directories(felitronics_reverb INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_reverb INTERFACE felitronics::core felitronics::eq)
target_compile_features(felitronics_reverb INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # Goldens rendered by the app's original stage, Freeverb's own properties, the stage contract, no allocation.
    add_executable(felitronics_reverb_tests tests/AmpReverbTests.cpp)
    target_link_libraries(felitronics_reverb_tests PRIVATE felitronics::reverb)
    target_include_directories(felitronics_reverb_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_reverb_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_reverb_tests COMMAND felitronics_reverb_tests)
endif()

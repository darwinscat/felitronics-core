# felitronics::multiband — the generic MultibandProcessor (split into LR4 bands → per-band processor →
# recombine) + the MultibandCompressor product. Header-only, JUCE-free. Deps: core (DelayLine) + eq
# (MultibandSplitter) + dynamics (Compressor).
add_library(felitronics_multiband INTERFACE)
add_library(felitronics::multiband ALIAS felitronics_multiband)
target_include_directories(felitronics_multiband INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_multiband INTERFACE felitronics::core felitronics::eq felitronics::dynamics)
target_compile_features(felitronics_multiband INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_multiband_tests tests/MultibandTests.cpp)
    target_link_libraries(felitronics_multiband_tests PRIVATE felitronics::multiband)
    target_include_directories(felitronics_multiband_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_multiband_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_multiband_tests COMMAND felitronics_multiband_tests)
endif()

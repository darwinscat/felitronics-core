# felitronics::lineareq — a linear-phase rendering of an eq bank (zero-phase FIR from the composite
# magnitude → click-free partitioned convolution). Header-only, JUCE-free. Deps: core (Fft) + eq
# (the magnitude grid) + convolution (ConvolutionEngine). Kept OUT of the base `eq` module so a consumer
# that only wants the matched biquads doesn't pull in the FFT convolver.
add_library(felitronics_lineareq INTERFACE)
add_library(felitronics::lineareq ALIAS felitronics_lineareq)
target_include_directories(felitronics_lineareq INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_lineareq INTERFACE felitronics::core felitronics::eq felitronics::convolution)
target_compile_features(felitronics_lineareq INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_lineareq_tests tests/LinearPhaseEqTests.cpp)
    target_link_libraries(felitronics_lineareq_tests PRIVATE felitronics::lineareq)
    target_include_directories(felitronics_lineareq_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_lineareq_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_lineareq_tests COMMAND felitronics_lineareq_tests)
endif()

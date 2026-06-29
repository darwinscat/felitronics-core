# felitronics::convolution — zero-latency partitioned FFT convolution on the core/Fft seam. Header-only
# for the spike (scalar reference FFT); becomes a compiled target when SIMD/real FFT backends + the
# offline IR resampler land. Depends on core (the FFT seam).
add_library(felitronics_convolution INTERFACE)
add_library(felitronics::convolution ALIAS felitronics_convolution)
target_include_directories(felitronics_convolution INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_convolution INTERFACE felitronics::core)
target_compile_features(felitronics_convolution INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_convolution_tests tests/ConvolutionTests.cpp)
    target_link_libraries(felitronics_convolution_tests PRIVATE felitronics::convolution)
    target_include_directories(felitronics_convolution_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_convolution_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_convolution_tests COMMAND felitronics_convolution_tests)

    add_executable(felitronics_convolution_resampler_tests tests/ResamplerTests.cpp)
    target_link_libraries(felitronics_convolution_resampler_tests PRIVATE felitronics::convolution)
    target_include_directories(felitronics_convolution_resampler_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_convolution_resampler_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_convolution_resampler_tests COMMAND felitronics_convolution_resampler_tests)

    add_executable(felitronics_convolution_engine_tests tests/ConvolutionEngineTests.cpp)
    target_link_libraries(felitronics_convolution_engine_tests PRIVATE felitronics::convolution)
    target_include_directories(felitronics_convolution_engine_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_convolution_engine_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_convolution_engine_tests COMMAND felitronics_convolution_engine_tests)
endif()

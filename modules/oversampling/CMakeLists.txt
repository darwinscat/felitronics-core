# felitronics::oversampling — polyphase windowed-sinc FIR up/down sampler (for true-peak detection +
# alias-free nonlinear processing). Header-only scalar reference; a SIMD backend can replace it later.
add_library(felitronics_oversampling INTERFACE)
add_library(felitronics::oversampling ALIAS felitronics_oversampling)
target_include_directories(felitronics_oversampling INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_oversampling INTERFACE felitronics::core)
target_compile_features(felitronics_oversampling INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_oversampling_tests tests/OversamplingTests.cpp)
    target_link_libraries(felitronics_oversampling_tests PRIVATE felitronics::oversampling)
    target_include_directories(felitronics_oversampling_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_oversampling_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_oversampling_tests COMMAND felitronics_oversampling_tests)
endif()

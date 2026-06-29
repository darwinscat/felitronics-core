# felitronics::limiter — brickwall true-peak limiter (oversample → limit → downsample; guaranteed
# true-peak ceiling). Header-only; depends on core + oversampling.
add_library(felitronics_limiter INTERFACE)
add_library(felitronics::limiter ALIAS felitronics_limiter)
target_include_directories(felitronics_limiter INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_limiter INTERFACE felitronics::core felitronics::oversampling)
target_compile_features(felitronics_limiter INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_limiter_tests tests/LimiterTests.cpp)
    target_link_libraries(felitronics_limiter_tests PRIVATE felitronics::limiter)
    target_include_directories(felitronics_limiter_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_limiter_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_limiter_tests COMMAND felitronics_limiter_tests)
endif()

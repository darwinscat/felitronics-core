# felitronics::delay — the amp plugins' tape-style echo (TapeDelay) and its sync/DARK parameter law.
# Header-only, JUCE-free, depends on felitronics::core only.
add_library(felitronics_delay INTERFACE)
add_library(felitronics::delay ALIAS felitronics_delay)
target_include_directories(felitronics_delay INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_delay INTERFACE felitronics::core)
target_compile_features(felitronics_delay INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # Goldens rendered by the app's original stage, the stage contract (laws 8/11a/11c), poison, no allocation.
    add_executable(felitronics_delay_tests tests/TapeDelayTests.cpp)
    target_link_libraries(felitronics_delay_tests PRIVATE felitronics::delay)
    target_include_directories(felitronics_delay_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_delay_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_delay_tests COMMAND felitronics_delay_tests)
endif()

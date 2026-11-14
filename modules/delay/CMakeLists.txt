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

    # The shipping class at its default head layout against a frozen copy of the one-head class, byte for byte.
    add_executable(felitronics_delay_parity_tests tests/TapeDelayParityTests.cpp)
    target_link_libraries(felitronics_delay_parity_tests PRIVATE felitronics::delay)
    target_include_directories(felitronics_delay_parity_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_delay_parity_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_delay_parity_tests COMMAND felitronics_delay_parity_tests)

    # The playback heads: placement, the figure, the pan law, the ramps and glides, re-slicing at event samples, law 8,
    # the tail, the setter's bounds, against oracles computed outside the class (the frozen one-head class among them).
    add_executable(felitronics_delay_heads_tests tests/TapeDelayHeadsTests.cpp)
    target_link_libraries(felitronics_delay_heads_tests PRIVATE felitronics::delay)
    target_include_directories(felitronics_delay_heads_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_delay_heads_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_delay_heads_tests COMMAND felitronics_delay_heads_tests)
endif()

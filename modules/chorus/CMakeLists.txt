# felitronics::chorus — the bass chorus (BassChorus) and its one-slider law.
# Header-only, JUCE-free, depends on felitronics::core and felitronics::eq (the crossover's DetMath Svf).
add_library(felitronics_chorus INTERFACE)
add_library(felitronics::chorus ALIAS felitronics_chorus)
target_include_directories(felitronics_chorus INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_chorus INTERFACE felitronics::core felitronics::eq)
target_compile_features(felitronics_chorus INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The stage contract (the law-11 properties, laws 8/8a/11a/11c), amount 0 = bit-exact dry, the split's
    # complement, mono lows, the measured pitch deviation, no click on any knob, poison, no allocation.
    add_executable(felitronics_chorus_tests tests/BassChorusTests.cpp)
    target_link_libraries(felitronics_chorus_tests PRIVATE felitronics::chorus)
    target_include_directories(felitronics_chorus_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_chorus_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_chorus_tests COMMAND felitronics_chorus_tests)

    # The listening tool: a bass DI (the WAVs given, or a synthesised line) through the stage at several
    # slider positions, written as 24-bit WAVs. Built so it cannot rot; run by hand.
    if(NOT EMSCRIPTEN)
        add_executable(felitronics_chorus_render tools/ChorusRender.cpp)
        target_link_libraries(felitronics_chorus_render PRIVATE felitronics::chorus felitronics::io)
        target_compile_features(felitronics_chorus_render PRIVATE cxx_std_20)
    endif()
endif()

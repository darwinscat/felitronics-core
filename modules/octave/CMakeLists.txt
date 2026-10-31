# felitronics::octave — the bass amp's sub-octave (Octaver) and its one-slider law.
# Header-only, JUCE-free; depends on core, eq (Svf) and dynamics (EnvelopeFollower).
add_library(felitronics_octave INTERFACE)
add_library(felitronics::octave ALIAS felitronics_octave)
target_include_directories(felitronics_octave INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_octave INTERFACE felitronics::core felitronics::eq felitronics::dynamics)
target_compile_features(felitronics_octave INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The tracker's own claims (pitch, lock, decay, harmonics, silence, latency), the stage contract (the law-11
    # properties, pause, reset, poison, law 8), the slider law, no allocation.
    add_executable(felitronics_octave_tests tests/OctaverTests.cpp)
    target_link_libraries(felitronics_octave_tests PRIVATE felitronics::octave felitronics::test_support)
    target_compile_features(felitronics_octave_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_octave_tests COMMAND felitronics_octave_tests)

    # Offline listening renders: a bass DI (or a synthesised bass line) through the slider law, 24-bit WAVs.
    add_executable(fcore_octave_render tools/octave_render.cpp)
    target_link_libraries(fcore_octave_render PRIVATE felitronics::octave felitronics::io)
    target_compile_features(fcore_octave_render PRIVATE cxx_std_20)
endif()

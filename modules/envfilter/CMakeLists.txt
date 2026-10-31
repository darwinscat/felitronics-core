# felitronics::envfilter — the bass envelope filter (EnvelopeFilter) and its one-slider law.
# Header-only, JUCE-free, depends on core, dynamics (the follower) and eq (the split's Svf).
add_library(felitronics_envfilter INTERFACE)
add_library(felitronics::envfilter ALIAS felitronics_envfilter)
target_include_directories(felitronics_envfilter INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_envfilter INTERFACE felitronics::core felitronics::dynamics felitronics::eq)
target_compile_features(felitronics_envfilter INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The sweep tracks the envelope, stability at the edges, no zipper, the clean lows, the slider law and the
    # stage contract (law 11 properties, 8a re-slicing, 11a/11c, poison, law 8, no allocation).
    add_executable(felitronics_envfilter_tests tests/EnvelopeFilterTests.cpp)
    target_link_libraries(felitronics_envfilter_tests PRIVATE felitronics::envfilter felitronics::test_support)
    target_compile_features(felitronics_envfilter_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_envfilter_tests COMMAND felitronics_envfilter_tests)
endif()

# A listening tool, not a test: renders a bass DI through the stage at several slider positions to 24-bit WAVs.
if(PROJECT_IS_TOP_LEVEL AND NOT EMSCRIPTEN)
    add_executable(felitronics_envfilter_render tools/EnvelopeFilterRender.cpp)
    target_link_libraries(felitronics_envfilter_render PRIVATE felitronics::envfilter felitronics::io)
    target_compile_features(felitronics_envfilter_render PRIVATE cxx_std_20)
endif()

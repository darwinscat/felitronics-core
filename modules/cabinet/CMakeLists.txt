# felitronics::cabinet — an amp's cabinet stage: one IR shot, what the player does to it (its tail trimmed, its cuts,
# its sign) baked into its taps, levelled to the guitar's reference or the bass's, resampled to the session's rate and
# handed to convolution::CabConvolver (zero latency, a click-free swap); the drain and the poison guard around it.
# Composition, so a module of its own: depends on core, eq (the baked cuts), measurement (the reference gains) and
# convolution. Header-only. Moved from OrbitAmp (its DSP half; decoding a file stays the caller's).
add_library(felitronics_cabinet INTERFACE)
add_library(felitronics::cabinet ALIAS felitronics_cabinet)
target_include_directories(felitronics_cabinet INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_cabinet INTERFACE felitronics::core felitronics::eq felitronics::measurement felitronics::convolution)
target_compile_features(felitronics_cabinet INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_cabinet_tests tests/CabinetIrTests.cpp)
    target_link_libraries(felitronics_cabinet_tests PRIVATE felitronics::cabinet)
    target_include_directories(felitronics_cabinet_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_cabinet_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_cabinet_tests COMMAND felitronics_cabinet_tests)
endif()

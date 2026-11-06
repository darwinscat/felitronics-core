# felitronics::basscomp — the bass amp's COMP: core's compressor and the bass's own characters (dynamics::CharacterComp)
# behind one blend, a fixed key high pass (eq::CutBiquad) and an automatic, level-neutral make-up; the standby is the
# house crossfade (core::BypassFade). Composition, so a module of its own (as deesser and envfilter are): dynamics
# itself stays independent of eq. Header-only; depends on core, eq and dynamics. Moved from OrbitAmp.
add_library(felitronics_basscomp INTERFACE)
add_library(felitronics::basscomp ALIAS felitronics_basscomp)
target_include_directories(felitronics_basscomp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_basscomp INTERFACE felitronics::core felitronics::eq felitronics::dynamics)
target_compile_features(felitronics_basscomp INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_basscomp_tests tests/BassCompTests.cpp)
    target_link_libraries(felitronics_basscomp_tests PRIVATE felitronics::basscomp)
    target_include_directories(felitronics_basscomp_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_basscomp_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_basscomp_tests COMMAND felitronics_basscomp_tests)
endif()

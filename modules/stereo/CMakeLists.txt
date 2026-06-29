# felitronics::stereo — stereo-field tools: the Mid/Side matrix + MonoBass (elliptical-EQ / bass mono-maker
# for vinyl + mastering). Header-only, JUCE-free. Depends on core + eq (the Cytomic SVF supplies the LR4
# low-pass + all-pass sections).
add_library(felitronics_stereo INTERFACE)
add_library(felitronics::stereo ALIAS felitronics_stereo)
target_include_directories(felitronics_stereo INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_stereo INTERFACE felitronics::core felitronics::eq)
target_compile_features(felitronics_stereo INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_stereo_tests tests/StereoTests.cpp)
    target_link_libraries(felitronics_stereo_tests PRIVATE felitronics::stereo)
    target_include_directories(felitronics_stereo_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_stereo_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_stereo_tests COMMAND felitronics_stereo_tests)
endif()

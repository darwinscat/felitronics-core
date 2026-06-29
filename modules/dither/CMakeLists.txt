# felitronics::dither — TPDF dither + error-feedback noise shaping, the last stage before bit-depth
# reduction on export (float → 16/20/24-bit). Header-only, JUCE-free. Depends only on core.
add_library(felitronics_dither INTERFACE)
add_library(felitronics::dither ALIAS felitronics_dither)
target_include_directories(felitronics_dither INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_dither INTERFACE felitronics::core)
target_compile_features(felitronics_dither INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_dither_tests tests/DitherTests.cpp)
    target_link_libraries(felitronics_dither_tests PRIVATE felitronics::dither)
    target_include_directories(felitronics_dither_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_dither_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_dither_tests COMMAND felitronics_dither_tests)
endif()

# felitronics::text — JUCE-free UTF-8 → ASCII romanization for building filesystem-safe slugs
# (German umlauts / ß + Russian Cyrillic → ASCII). Header-only INTERFACE; zero deps (not even
# felitronics::core). MESSAGE-THREAD-ONLY (allocates). Extracted from OrbitCapture NAM's Matrix.h so
# every Darwin's Cat product shares ONE romanization table.
add_library(felitronics_text INTERFACE)
add_library(felitronics::text ALIAS felitronics_text)
target_include_directories(felitronics_text INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(felitronics_text INTERFACE cxx_std_17)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_text_tests tests/TranslitTests.cpp)
    target_link_libraries(felitronics_text_tests PRIVATE felitronics::text)
    target_include_directories(felitronics_text_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_text_tests PRIVATE cxx_std_17)
    add_test(NAME felitronics_text_tests COMMAND felitronics_text_tests)
endif()

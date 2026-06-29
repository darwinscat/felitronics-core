# felitronics::deesser — a de-esser with two topologies (surgical dynamic-EQ + split-band duck). Header-only,
# JUCE-free. Composes eq (Crossover2/Svf) + dynamics (detector/gain) + dynamiceq (the DynamicEq mode).
add_library(felitronics_deesser INTERFACE)
add_library(felitronics::deesser ALIAS felitronics_deesser)
target_include_directories(felitronics_deesser INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_deesser INTERFACE felitronics::core felitronics::eq felitronics::dynamics felitronics::dynamiceq)
target_compile_features(felitronics_deesser INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_deesser_tests tests/DeEsserTests.cpp)
    target_link_libraries(felitronics_deesser_tests PRIVATE felitronics::deesser)
    target_include_directories(felitronics_deesser_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_deesser_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_deesser_tests COMMAND felitronics_deesser_tests)
endif()

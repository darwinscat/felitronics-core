# felitronics::dynamiceq — a Dynamic EQ band (eq::Svf gain driven by a dynamics detector). Header-only,
# JUCE-free. Deps: core + eq (the filters) + dynamics (detector / gain computer / ballistics). It composes
# eq + dynamics, which is why it isn't in either of them.
add_library(felitronics_dynamiceq INTERFACE)
add_library(felitronics::dynamiceq ALIAS felitronics_dynamiceq)
target_include_directories(felitronics_dynamiceq INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_dynamiceq INTERFACE felitronics::core felitronics::eq felitronics::dynamics)
target_compile_features(felitronics_dynamiceq INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_dynamiceq_tests tests/DynamicEqTests.cpp)
    target_link_libraries(felitronics_dynamiceq_tests PRIVATE felitronics::dynamiceq)
    target_include_directories(felitronics_dynamiceq_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_dynamiceq_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_dynamiceq_tests COMMAND felitronics_dynamiceq_tests)
endif()

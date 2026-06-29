# felitronics::neural — the process-only inference seam (Inference concept) + NeuralStage swap-safe
# holder. Header-only, JUCE-free AND NAM-free (the NAM/Eigen backend + model loading live in the
# adapter). Depends only on core.
add_library(felitronics_neural INTERFACE)
add_library(felitronics::neural ALIAS felitronics_neural)
target_include_directories(felitronics_neural INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_neural INTERFACE felitronics::core)
target_compile_features(felitronics_neural INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_neural_tests tests/NeuralTests.cpp)
    target_link_libraries(felitronics_neural_tests PRIVATE felitronics::neural)
    target_include_directories(felitronics_neural_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_neural_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_neural_tests COMMAND felitronics_neural_tests)
endif()

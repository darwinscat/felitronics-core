# felitronics::measurement — OFFLINE IR-measurement math: ESS sweep + Farina inverse, Farina
# deconvolution (latency-absorbing onset past the harmonic region), IR post (onset/trim/normalize),
# a pre-deconv capture-quality gate, and multi-mic common-onset alignment. Extracted from OrbitCapture's
# capture engine (written portable by design). Header-only INTERFACE; depends on core. MESSAGE-THREAD-
# ONLY (double, allocates) — this is NOT the RT FFT seam (RT consumers use felitronics::core::fft).
add_library(felitronics_measurement INTERFACE)
add_library(felitronics::measurement ALIAS felitronics_measurement)
target_include_directories(felitronics_measurement INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_measurement INTERFACE felitronics::core)
target_compile_features(felitronics_measurement INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_measurement_tests tests/MeasurementTests.cpp)
    target_link_libraries(felitronics_measurement_tests PRIVATE felitronics::measurement)
    target_include_directories(felitronics_measurement_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_measurement_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_measurement_tests COMMAND felitronics_measurement_tests)

    add_executable(felitronics_measurement_adversarial_tests tests/AdversarialTests.cpp)
    target_link_libraries(felitronics_measurement_adversarial_tests PRIVATE felitronics::measurement)
    target_include_directories(felitronics_measurement_adversarial_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_measurement_adversarial_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_measurement_adversarial_tests COMMAND felitronics_measurement_adversarial_tests)

    add_executable(felitronics_modelguess_tests tests/ModelGuessTests.cpp)
    target_link_libraries(felitronics_modelguess_tests PRIVATE felitronics::measurement)
    target_include_directories(felitronics_modelguess_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_modelguess_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_modelguess_tests COMMAND felitronics_modelguess_tests)
endif()

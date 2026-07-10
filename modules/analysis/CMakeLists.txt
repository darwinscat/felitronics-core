# felitronics::analysis — SpectrumTap (lock-free SPSC analyser tap) today; correlation / LUFS and the
# FFT-via-seam consumers land later. Header-only; depends on core.
add_library(felitronics_analysis INTERFACE)
add_library(felitronics::analysis ALIAS felitronics_analysis)
target_include_directories(felitronics_analysis INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_analysis INTERFACE felitronics::core)
target_compile_features(felitronics_analysis INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_analysis_tests tests/AnalysisTests.cpp)
    target_link_libraries(felitronics_analysis_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_analysis_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_analysis_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_analysis_tests COMMAND felitronics_analysis_tests)

    add_executable(felitronics_truepeak_tests tests/TruePeakMeterTests.cpp)
    target_link_libraries(felitronics_truepeak_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_truepeak_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_truepeak_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_truepeak_tests COMMAND felitronics_truepeak_tests)

    add_executable(felitronics_spectrum_curve_tests tests/SpectrumCurveTests.cpp)
    target_link_libraries(felitronics_spectrum_curve_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_spectrum_curve_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_spectrum_curve_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_spectrum_curve_tests COMMAND felitronics_spectrum_curve_tests)
endif()

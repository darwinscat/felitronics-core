# felitronics::analysis — RT metering, the analyser taps and the spectrum panes (SpectrumPane, MultiResSpectrumPane);
# offline display curves under ::offline. Header-only; depends on core.
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

    add_executable(felitronics_loudness_conformance_tests tests/LoudnessConformanceTests.cpp)   # EBU Tech 3341 Table 1 + gates + K + chunk invariance
    target_link_libraries(felitronics_loudness_conformance_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_loudness_conformance_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_loudness_conformance_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_loudness_conformance_tests COMMAND felitronics_loudness_conformance_tests)

    add_executable(felitronics_truepeak_conformance_tests tests/TruePeakConformanceTests.cpp)   # EBU Tech 3341 §2.6 true-peak tests 15-19
    target_link_libraries(felitronics_truepeak_conformance_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_truepeak_conformance_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_truepeak_conformance_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_truepeak_conformance_tests COMMAND felitronics_truepeak_conformance_tests)

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

    add_executable(felitronics_spectrum_pane_perf_tests tests/SpectrumPanePerfTests.cpp)
    target_link_libraries(felitronics_spectrum_pane_perf_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_spectrum_pane_perf_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_spectrum_pane_perf_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_spectrum_pane_perf_tests COMMAND felitronics_spectrum_pane_perf_tests)

    add_executable(felitronics_multires_fast_pane_tests tests/MultiResFastPaneTests.cpp)   # paired NULL: the fast sibling vs the pane it copies
    target_link_libraries(felitronics_multires_fast_pane_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_multires_fast_pane_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_multires_fast_pane_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_multires_fast_pane_tests COMMAND felitronics_multires_fast_pane_tests)

    # WaveformPeaks + StereoColumns — ports of the site's audio-peaks.js / stereo-meter.js, pinned to the spec's own
    # numbers. Built with the tree's own `on`, which is where the natural removal of the stereo products' volatile pin
    # (`mid += m*m`) fuses on arm64 and fails witnesses W5/W6f — see the test file's header for the measured matrix.
    add_executable(felitronics_waveform_shape_tests tests/WaveformShapeTests.cpp)
    target_link_libraries(felitronics_waveform_shape_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_waveform_shape_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_waveform_shape_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_waveform_shape_tests COMMAND felitronics_waveform_shape_tests)

    add_executable(felitronics_multires_spectrum_tests tests/MultiResSpectrumTests.cpp)
    target_link_libraries(felitronics_multires_spectrum_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_multires_spectrum_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_multires_spectrum_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_multires_spectrum_tests COMMAND felitronics_multires_spectrum_tests)
endif()

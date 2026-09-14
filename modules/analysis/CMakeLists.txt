# felitronics::analysis — RT metering, the analyser taps and the spectrum panes (SpectrumPane, MultiResSpectrumPane);
# offline display curves under ::offline. Header-only; depends on core, and on oversampling for the REFERENCE
# true-peak meter (ReferenceTruePeakMeter), which is the reference because it is that oversampler.
add_library(felitronics_analysis INTERFACE)
add_library(felitronics::analysis ALIAS felitronics_analysis)
target_include_directories(felitronics_analysis INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_analysis INTERFACE felitronics::core felitronics::oversampling)
target_compile_features(felitronics_analysis INTERFACE cxx_std_20)

# felitronics::analysis_offline — the link surface for the OFFLINE analyzers (ProgrammeReport,
# SourceForensics, HumDetector, LowEnd, BandBursts), which compose over the eq/stereo primitives
# (eq::Crossover2 is the LR4, eq::Svf the biquad, stereo::MidSide the M/S pair) instead of hand-rolling
# a fifth private copy of each. It CANNOT be a dependency of `analysis` itself: `eq` links
# `felitronics::analysis` (modules/eq/CMakeLists.txt:7 — EqTypes.h:18 includes analysis/SpectrumTap.h)
# and `stereo` links `eq`, so analysis -> eq would close a CYCLE. Separating the link surface keeps the
# graph acyclic; the headers stay under analysis/ and the namespace stays felitronics::analysis.
# Declared BEFORE the test block below so the analyzers' test targets can link it, and forward-referencing
# felitronics::eq / ::stereo (added later in the top-level CMakeLists) is resolved at generate time.
add_library(felitronics_analysis_offline INTERFACE)
add_library(felitronics::analysis_offline ALIAS felitronics_analysis_offline)
target_link_libraries(felitronics_analysis_offline INTERFACE
    felitronics::analysis felitronics::eq felitronics::stereo)
target_compile_features(felitronics_analysis_offline INTERFACE cxx_std_20)

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

    add_executable(felitronics_spectrum_frames_tests tests/SpectrumFramesTests.cpp)   # the shared frame producer: DFT null + law 8a
    target_link_libraries(felitronics_spectrum_frames_tests PRIVATE felitronics::analysis_offline)
    target_include_directories(felitronics_spectrum_frames_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_spectrum_frames_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_spectrum_frames_tests COMMAND felitronics_spectrum_frames_tests)

    add_executable(felitronics_clip_detector_tests tests/ClipDetectorTests.cpp)   # clipping by construction + a whole-file reference null
    target_link_libraries(felitronics_clip_detector_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_clip_detector_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_clip_detector_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_clip_detector_tests COMMAND felitronics_clip_detector_tests)

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

    add_executable(felitronics_reference_truepeak_tests tests/ReferenceTruePeakMeterTests.cpp)   # bit-for-bit the 4x/32 reference
    target_link_libraries(felitronics_reference_truepeak_tests PRIVATE felitronics::analysis)
    target_include_directories(felitronics_reference_truepeak_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_reference_truepeak_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_reference_truepeak_tests COMMAND felitronics_reference_truepeak_tests)

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

    # ProgrammeReport (P72) — the OFFLINE whole-programme report. Links felitronics::analysis_offline, not
    # felitronics::analysis: it composes over eq::Crossover2 (the LR4) and the stereo primitives, and adding
    # those to `analysis` itself would close a cycle (see the comment at that target above). The suite's own
    # whole-file reference drives eq::Crossover2 and analysis::KWeightingFilter directly, so it needs the
    # same surface.
    add_executable(felitronics_programme_report_tests tests/ProgrammeReportTests.cpp)
    target_link_libraries(felitronics_programme_report_tests PRIVATE felitronics::analysis_offline)
    target_include_directories(felitronics_programme_report_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_programme_report_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_programme_report_tests COMMAND felitronics_programme_report_tests)
    # HumDetector — mains hum against the bass notes that look like it (49.0 / 58.3 / 61.7 Hz), plus law 8a
    # over the whole report. Links analysis_offline: the offline analyzers' link surface (see the note above).
    add_executable(felitronics_hum_detector_tests tests/HumDetectorTests.cpp)
    target_link_libraries(felitronics_hum_detector_tests PRIVATE felitronics::analysis_offline)
    target_include_directories(felitronics_hum_detector_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_hum_detector_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_hum_detector_tests COMMAND felitronics_hum_detector_tests)
    # LowEnd — wide bass (LR4 + M/S over the 10 ms grid) and the dominant low note (semitone fold of the
    # shared frame producer). Links analysis_offline: it composes over eq::Crossover2 and stereo::MidSide,
    # which `analysis` itself cannot depend on without closing a cycle (see the target's own note above).
    add_executable(felitronics_low_end_tests tests/LowEndTests.cpp)
    target_link_libraries(felitronics_low_end_tests PRIVATE felitronics::analysis_offline)
    target_include_directories(felitronics_low_end_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_low_end_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_low_end_tests COMMAND felitronics_low_end_tests)
    # BandBursts (P76) — an offline burst detector: energy in a band against the MEDIAN of a long trailing
    # ring of hops. Links analysis_offline for eq::Crossover2 (the LR4 that makes the band).
    add_executable(felitronics_band_bursts_tests tests/BandBurstsTests.cpp)
    target_link_libraries(felitronics_band_bursts_tests PRIVATE felitronics::analysis_offline)
    target_include_directories(felitronics_band_bursts_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_band_bursts_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_band_bursts_tests COMMAND felitronics_band_bursts_tests)
endif()

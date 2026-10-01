# felitronics::codecgrid — OFFLINE measurement: the frame grid a lossy codec leaves in decoded PCM. A stretch of
# audio is analysed again with a codec's own transform (the MPEG-1 Layer III hybrid filterbank, the AAC and CELT
# MDCTs) at every frame offset it could have been coded at; at the encoder's offset the quantiser's zeros come
# back, one sample off they do not. Header-only INTERFACE; depends on core. Not a codec: nothing here encodes or
# decodes, and no bitstream is read.
add_library(felitronics_codecgrid INTERFACE)
add_library(felitronics::codecgrid ALIAS felitronics_codecgrid)
target_include_directories(felitronics_codecgrid INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_codecgrid INTERFACE felitronics::core)
target_compile_features(felitronics_codecgrid INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The MP3 analysis filterbank against the standard's formulas written out directly.
    add_executable(felitronics_codecgrid_mp3hybrid_tests tests/Mp3HybridTests.cpp)
    target_link_libraries(felitronics_codecgrid_mp3hybrid_tests PRIVATE felitronics::codecgrid)
    target_include_directories(felitronics_codecgrid_mp3hybrid_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_codecgrid_mp3hybrid_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_codecgrid_mp3hybrid_tests COMMAND felitronics_codecgrid_mp3hybrid_tests)

    # The curve: the level accumulator against a logarithm per coefficient, and the reading of a curve.
    add_executable(felitronics_codecgrid_curve_tests tests/GridCurveTests.cpp)
    target_link_libraries(felitronics_codecgrid_curve_tests PRIVATE felitronics::codecgrid)
    target_include_directories(felitronics_codecgrid_curve_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_codecgrid_curve_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_codecgrid_curve_tests COMMAND felitronics_codecgrid_curve_tests)

    # The complex FFT, the windows and the MDCT against their definitions.
    add_executable(felitronics_codecgrid_mdct_tests tests/MdctTests.cpp)
    target_link_libraries(felitronics_codecgrid_mdct_tests PRIVATE felitronics::codecgrid)
    target_include_directories(felitronics_codecgrid_mdct_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_codecgrid_mdct_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_codecgrid_mdct_tests COMMAND felitronics_codecgrid_mdct_tests)

    # The scan on programmes made to carry a grid (a transform, zeros, the inverse — no codec, no audio file),
    # and the back resampler.
    add_executable(felitronics_codecgrid_scan_tests tests/GridScanTests.cpp)
    target_link_libraries(felitronics_codecgrid_scan_tests PRIVATE felitronics::codecgrid)
    target_include_directories(felitronics_codecgrid_scan_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_codecgrid_scan_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_codecgrid_scan_tests COMMAND felitronics_codecgrid_scan_tests)

    # The detector's verdicts on programmes made to have one property each.
    add_executable(felitronics_codecgrid_detector_tests tests/CodecGridDetectorTests.cpp)
    target_link_libraries(felitronics_codecgrid_detector_tests PRIVATE felitronics::codecgrid)
    target_include_directories(felitronics_codecgrid_detector_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_codecgrid_detector_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_codecgrid_detector_tests COMMAND felitronics_codecgrid_detector_tests)

    # The detector's form: the call contract (law 11), split invariance, the published storage, no allocation.
    add_executable(felitronics_codecgrid_law_tests tests/CodecGridLawTests.cpp)
    target_link_libraries(felitronics_codecgrid_law_tests PRIVATE felitronics::codecgrid)
    target_include_directories(felitronics_codecgrid_law_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_codecgrid_law_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_codecgrid_law_tests COMMAND felitronics_codecgrid_law_tests)
endif()

# felitronics::poweramp — the white-box tube power-amp trunk, lifted verbatim from OrbitCab:
# SagEnvelope (feedforward dual-TC supply sag) + TubeStage (PP/SE composite tube transfer) +
# PowerAmpStage (the oversampled chain: virtual load → drive/sag rail → waveshaper + DC-block +
# output transformer in the OS domain → MID bell → NFB shelves → drive-comp/output). Header-only;
# the PRODUCT supplies the per-tube voicing data (Voicing) + control values (Params) — OrbitCab's
# kTubeVoicings preset table and its router/param mapping stay in the product.
add_library(felitronics_poweramp INTERFACE)
add_library(felitronics::poweramp ALIAS felitronics_poweramp)
target_include_directories(felitronics_poweramp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_poweramp INTERFACE
    felitronics::core felitronics::eq felitronics::oversampling felitronics::saturation)
target_compile_features(felitronics_poweramp INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    # The MANIACAL golden battery ported from OrbitCab's tests/PowerAmpCoreGolden.cpp — it pins the
    # stage's numeric fingerprint (same thresholds, same hostile schedules, same voicing constants).
    add_executable(felitronics_poweramp_golden_tests tests/PowerAmpGoldenTests.cpp)
    target_link_libraries(felitronics_poweramp_golden_tests PRIVATE felitronics::poweramp)
    target_include_directories(felitronics_poweramp_golden_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_poweramp_golden_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_poweramp_golden_tests COMMAND felitronics_poweramp_golden_tests)

    # THEORY-FIRST falsification suite — derives every bound from amplifier/DSP theory + the headers'
    # documented contracts (never tuned to pass) and tries to BREAK the stage. Documented deviations
    # (the ~1e-2 "B9" feel-layer block-size discrepancy) are PINNED here, not fixed.
    add_executable(felitronics_poweramp_theory_tests tests/PowerAmpTheoryTests.cpp)
    target_link_libraries(felitronics_poweramp_theory_tests PRIVATE felitronics::poweramp)
    target_include_directories(felitronics_poweramp_theory_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_poweramp_theory_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_poweramp_theory_tests COMMAND felitronics_poweramp_theory_tests)
endif()

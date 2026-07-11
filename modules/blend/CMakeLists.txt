# felitronics::blend — the mic-blend engine: mixer strip / Master parameter model, pure DSP kernels
# (biquad/one-pole/Butterworth-slope cascade/Hilbert phase rotation/fractional shift), and the blend
# engine itself (solo/mute, per-strip filters+shift, Σ(gain·strip), Master bus). Promoted VERBATIM from
# OrbitCapture's app-local mixer (model/MixModel.h + core/BlendKernels.h + core/BlendEngine.h).
# Header-only INTERFACE; OFFLINE (double coefficient design, float sample buffers) — depends on core
# (Math::kPi, OfflineFft) + analysis.
add_library(felitronics_blend INTERFACE)
add_library(felitronics::blend ALIAS felitronics_blend)
target_include_directories(felitronics_blend INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_blend INTERFACE felitronics::core felitronics::analysis)
target_compile_features(felitronics_blend INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_blend_tests tests/BlendTests.cpp)
    target_link_libraries(felitronics_blend_tests PRIVATE felitronics::blend)
    target_include_directories(felitronics_blend_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_blend_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_blend_tests COMMAND felitronics_blend_tests)

    add_executable(felitronics_blend_overlay_tests tests/OverlayTests.cpp)
    target_link_libraries(felitronics_blend_overlay_tests PRIVATE felitronics::blend)
    target_include_directories(felitronics_blend_overlay_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_blend_overlay_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_blend_overlay_tests COMMAND felitronics_blend_overlay_tests)
endif()

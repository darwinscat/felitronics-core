# felitronics-core dev tools (built only when this is the top-level project). Not shipped; used to
# validate the DSP against reference implementations (ffmpeg) on real audio.
# The FP-contraction contract for everything that goes through fcore_probe.h. Baseline wasm has NO scalar FMA
# instruction, so emscripten cannot contract a*b+c; a contracting native build therefore disagrees with wasm —
# measured at 6.75e-7 dB on the true-peak answer, which is 675x the tolerance the P0 spike is held to. GCC is
# the sharper edge here: it defaults to -ffp-contract=fast in gnu++ modes, so a Linux/GCC reference contracts
# even where clang would not. felitronics_probe_tests asserts the flag actually took effect, so dropping this
# line turns a silent parity failure into a red test.
# (Not applied on MSVC: a Windows build is not part of the parity contract, and /fp: has no equivalent knob.)
set(FCORE_NO_FP_CONTRACT $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-ffp-contract=off>)

add_executable(fcore_measure fcore_measure.cpp)
target_link_libraries(fcore_measure PRIVATE felitronics::analysis felitronics::oversampling)
target_compile_features(fcore_measure PRIVATE cxx_std_20)
target_compile_options(fcore_measure PRIVATE ${FCORE_NO_FP_CONTRACT})

if(FELITRONICS_BUILD_TESTS)
    # fcore::Probe is the body shared verbatim with the wasm shim, so it is the one piece of tools/ that
    # carries a correctness claim rather than convenience — and it gets a real suite.
    add_executable(felitronics_probe_tests tests/ProbeTests.cpp)
    target_link_libraries(felitronics_probe_tests PRIVATE felitronics::analysis felitronics::oversampling)
    target_include_directories(felitronics_probe_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_probe_tests PRIVATE cxx_std_20)
    target_compile_options(felitronics_probe_tests PRIVATE ${FCORE_NO_FP_CONTRACT})
    add_test(NAME felitronics_probe_tests COMMAND felitronics_probe_tests)

    # The wasm C ABI is the one surface fed by untrusted input (whatever a page's JavaScript computed), so its
    # guards are tested — natively, because fc_probe.cpp compiles anywhere (EMSCRIPTEN_KEEPALIVE degrades to a
    # plain extern "C"). That puts the validation and the planar addressing under ctest, ASan and UBSan.
    add_executable(felitronics_abi_tests tests/AbiTests.cpp wasm/fc_probe.cpp)
    target_link_libraries(felitronics_abi_tests PRIVATE felitronics::analysis felitronics::oversampling)
    target_include_directories(felitronics_abi_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_abi_tests PRIVATE cxx_std_20)
    target_compile_options(felitronics_abi_tests PRIVATE ${FCORE_NO_FP_CONTRACT})
    add_test(NAME felitronics_abi_tests COMMAND felitronics_abi_tests)
endif()

add_executable(fcore_convolve fcore_convolve.cpp)
target_link_libraries(fcore_convolve PRIVATE felitronics::convolution)
target_compile_features(fcore_convolve PRIVATE cxx_std_20)

# fftbench — the #1-debt acceptance benchmark (convolver %RT: old P>=block vs new fixed P=128, scalar vs
# pffft). Build Release; the pffft column appears only with -DFELITRONICS_WITH_PFFFT=ON.
add_executable(fcore_fftbench fftbench.cpp)
target_link_libraries(fcore_fftbench PRIVATE felitronics::convolution)
target_compile_features(fcore_fftbench PRIVATE cxx_std_20)
if(FELITRONICS_WITH_PFFFT)
    target_link_libraries(fcore_fftbench PRIVATE felitronics::fftpffft)
    target_compile_definitions(fcore_fftbench PRIVATE FELITRONICS_WITH_PFFFT)
endif()

# Optional JUCE/vDSP FFT reference column — FetchContent's JUCE INTO THE BENCH ONLY (never the core/product).
option(FELITRONICS_BENCH_JUCE "Fetch JUCE and add a JUCE-FFT reference column to fcore_fftbench" OFF)
if(FELITRONICS_BENCH_JUCE)
    include(FetchContent)
    FetchContent_Declare(juce GIT_REPOSITORY https://github.com/juce-framework/JUCE.git GIT_TAG 8.0.4 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(juce)
    target_link_libraries(fcore_fftbench PRIVATE juce::juce_dsp juce::juce_events)
    target_compile_definitions(fcore_fftbench PRIVATE
        FELITRONICS_BENCH_JUCE=1 JUCE_USE_CURL=0 JUCE_WEB_BROWSER=0 JUCE_STANDALONE_APPLICATION=1)
endif()

# eqviz — renders PNG "quality pictures" (heatmaps + curve sheets) of every eq filter type against
# independent analog references. Offline; catches parameter-corner bugs eyes-first (see tools/eqviz).
add_executable(eqviz eqviz/eqviz.cpp)
target_link_libraries(eqviz PRIVATE felitronics::eq)
target_compile_features(eqviz PRIVATE cxx_std_20)

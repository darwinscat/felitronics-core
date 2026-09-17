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
target_link_libraries(fcore_measure PRIVATE felitronics::analysis_offline felitronics::oversampling)
# The OFFLINE analyzer modes (lowend, and the siblings landing beside it) compose over eq::Crossover2 /
# stereo::MidSide, which felitronics::analysis cannot itself depend on without closing a cycle — see the
# note on the target in modules/analysis/CMakeLists.txt. Added rather than substituted, so two branches
# that both need it resolve to one line instead of a conflict.
target_compile_features(fcore_measure PRIVATE cxx_std_20)
target_compile_options(fcore_measure PRIVATE ${FCORE_NO_FP_CONTRACT})

if(FELITRONICS_BUILD_TESTS)
    # fcore::Probe is a body shared verbatim with the wasm shim, so it is one of the two pieces of tools/ that
    # carry a correctness claim rather than convenience — and it gets a real suite. (The other is
    # fcore::ClipProbe, whose suite is felitronics_clips_exposure_tests, further down.)
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
    target_link_libraries(felitronics_abi_tests PRIVATE felitronics::analysis_offline felitronics::oversampling)
    target_include_directories(felitronics_abi_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_abi_tests PRIVATE cxx_std_20)
    target_compile_options(felitronics_abi_tests PRIVATE ${FCORE_NO_FP_CONTRACT})
    add_test(NAME felitronics_abi_tests COMMAND felitronics_abi_tests)

    # The EXPOSURE of analysis::ClipDetector (P71): fcore::ClipProbe, the `clips` text format and the clips
    # entry points of the wasm ABI. Same arrangement as felitronics_abi_tests one block up — fc_probe.cpp is
    # compiled natively, so the ABI's validation and packing run under ctest, ASan and UBSan.
    #
    # FCORE_NO_FP_CONTRACT IS LOAD-BEARING ON THIS TARGET, not inherited decoration. The repository declares
    # -ffp-contract=on (top-level CMakeLists, law 10), and the one expression of this measurement a fused
    # multiply-add can reach is ClipDetector's smooth-crest bound V = (tau + q)*stretch^2 + q — which gates a
    # strict >=, so a 1-ulp difference creates or destroys a whole RUN rather than moving a last digit. The
    # numbers this suite pins are the numbers `fcore_measure clips` and the wasm module produce, and both of
    # those are built contract-off; a suite built contract-on would be pinning a third answer. -fno-fast-math
    # is stated for the same reason it is stated in tools/wasm/build.sh: nothing here may be reassociated.
    # The P77 ABI surface under adversarial arguments — natively, so ASan/UBSan see it. fc_probe.cpp
    # compiles anywhere; EMSCRIPTEN_KEEPALIVE degrades to nothing off the wasm tier.
    add_executable(felitronics_analysis_abi_tests tests/AnalysisAbiTests.cpp wasm/fc_probe.cpp)
    target_link_libraries(felitronics_analysis_abi_tests PRIVATE felitronics::analysis_offline felitronics::oversampling)
    target_include_directories(felitronics_analysis_abi_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_features(felitronics_analysis_abi_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_analysis_abi_tests COMMAND felitronics_analysis_abi_tests)

    add_executable(felitronics_clips_exposure_tests tests/ClipsExposureTests.cpp wasm/fc_probe.cpp)
    target_link_libraries(felitronics_clips_exposure_tests PRIVATE felitronics::analysis_offline felitronics::oversampling)
    target_include_directories(felitronics_clips_exposure_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_clips_exposure_tests PRIVATE cxx_std_20)
    target_compile_options(felitronics_clips_exposure_tests PRIVATE ${FCORE_NO_FP_CONTRACT}
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-fno-fast-math>)
    add_test(NAME felitronics_clips_exposure_tests COMMAND felitronics_clips_exposure_tests)

    # The streaming surface of the same ABI (P120): the fc_stream_* handles over fcore::StreamProbe, fed in pieces,
    # against the same two instruments handed the whole buffer. Contract-off for the clips suite's reason — the
    # detector's smooth-crest bound decides whether a run exists, and both parity roads are built this way.
    add_executable(felitronics_stream_abi_tests tests/StreamAbiTests.cpp wasm/fc_probe.cpp)
    target_link_libraries(felitronics_stream_abi_tests PRIVATE felitronics::analysis_offline felitronics::oversampling)
    target_include_directories(felitronics_stream_abi_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_stream_abi_tests PRIVATE cxx_std_20)
    target_compile_options(felitronics_stream_abi_tests PRIVATE ${FCORE_NO_FP_CONTRACT}
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-fno-fast-math>)
    add_test(NAME felitronics_stream_abi_tests COMMAND felitronics_stream_abi_tests)

    # The mastering C ABI. Same arrangement as felitronics_abi_tests one block up — the shim compiles
    # natively, so its validation, ownership and mapping layer runs under ctest, ASan and UBSan. NOT
    # given FCORE_NO_FP_CONTRACT, for the reason spelled out at fcore_master above: this suite and the
    # `mastering` module suite must be the same arithmetic or a comparison between them means nothing.
    add_executable(felitronics_master_abi_tests tests/MasterAbiTests.cpp wasm/fc_master.cpp)
    target_link_libraries(felitronics_master_abi_tests PRIVATE felitronics::mastering felitronics::analysis)
    target_include_directories(felitronics_master_abi_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_master_abi_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_master_abi_tests COMMAND felitronics_master_abi_tests)

    # Re-entry poisons the module, for good — a binary of its own because the poison is permanent by contract and
    # the suite above spends its one poisoning on a call that never returned. It needs no exceptions, so the wasm
    # tier builds and runs it too.
    add_executable(felitronics_master_abi_reentry_tests tests/MasterAbiReentryTests.cpp wasm/fc_master.cpp)
    target_link_libraries(felitronics_master_abi_reentry_tests PRIVATE felitronics::mastering felitronics::analysis)
    target_include_directories(felitronics_master_abi_reentry_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_master_abi_reentry_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_master_abi_reentry_tests COMMAND felitronics_master_abi_reentry_tests)

    add_executable(felitronics_master_abi_progress_tests tests/MasterAbiProgressTests.cpp wasm/fc_master.cpp)
    target_link_libraries(felitronics_master_abi_progress_tests PRIVATE felitronics::mastering felitronics::analysis)
    target_include_directories(felitronics_master_abi_progress_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_master_abi_progress_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_master_abi_progress_tests COMMAND felitronics_master_abi_progress_tests)

    # The bit-exactness acceptance itself, run as a test rather than by hand: the same programme through
    # the C ABI and through a direct C++ call, one binary, one machine, compared bit for bit — plus
    # block-independence through the ABI and the reset discipline.
    add_test(NAME felitronics_master_abi_selftest COMMAND fcore_master selftest)

    # The JavaScript layout against the compiler, field by field (tools/wasm/layout-check.mjs). On the wasm tier the
    # argument is the emscripten build of fcore_master, so the offsets checked are wasm32's — the ones a page uses.
    # Needs node; a desktop row without it skips the test and says so rather than passing it.
    find_program(FELITRONICS_NODE node)
    if(FELITRONICS_NODE)
        add_test(NAME felitronics_master_layout_check
                 COMMAND ${FELITRONICS_NODE} ${CMAKE_CURRENT_SOURCE_DIR}/wasm/layout-check.mjs $<TARGET_FILE:fcore_master>)
        # And the two refusals the layout file makes for a page — the version gate and the write past a stamp.
        add_test(NAME felitronics_master_layout_gate_tests
                 COMMAND ${FELITRONICS_NODE} ${CMAKE_CURRENT_SOURCE_DIR}/wasm/layout-gate-test.mjs)
    else()
        message(STATUS "felitronics_master_layout_check: node not found — the JS layout is not checked on this row")
    endif()
endif()

# fcore_master — the C-ABI facade over `felitronics::mastering` (tools/fc_master_abi.h +
# tools/wasm/fc_master.cpp) and the native CLI that drives it. Two things to notice, both load-bearing:
#
#  1. IT DELIBERATELY DOES **NOT** TAKE `FCORE_NO_FP_CONTRACT`. Every other target in this file carries
#     `-ffp-contract=off` because it is part of the P0 probe's cross-toolchain contract; this one is the
#     reference for the C++ API the DESKTOP app links, so it takes the repository's declared
#     `-ffp-contract=on` (law 10 — the flag is STATED, never inherited by accident). Measured on one
#     source, the two flavours differ in 203269 of 288000 samples (worst 2.03e-6). Mixing them WITHIN a
#     binary is worse than either: two TUs built with different flags measured as AGREEING, because the
#     linker merged the header-only instantiations and whichever TU came first won — agreement by link
#     order, which is exactly the false-green `selftest` exists to prevent. So the whole target, both
#     halves of the comparison included, is built one way and that way is stated here.
#  2. THE CLI AND ITS ABI ARE SEPARATE TUs ON PURPOSE. `selftest` compares a render through the C entry
#     points against a direct C++ call; inlining the facade into the caller would make the "ABI" half
#     stop being an ABI at all.
add_executable(fcore_master fcore_master.cpp wasm/fc_master.cpp)
target_link_libraries(fcore_master PRIVATE felitronics::mastering felitronics::analysis)
target_include_directories(fcore_master PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(fcore_master PRIVATE cxx_std_20)

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

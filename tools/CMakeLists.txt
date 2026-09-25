# felitronics-core dev tools (built only when this is the top-level project). Not shipped: a convolution CLI, the
# convolver benchmark and the EQ quality pictures. The C ABIs of the mastering chain and the offline analyzers,
# their CLIs (fcore_master, fcore_measure) and the wasm build live in felitronics-mastering-core.

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

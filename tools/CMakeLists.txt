# felitronics-core dev tools (built only when this is the top-level project). Not shipped; used to
# validate the DSP against reference implementations (ffmpeg) on real audio.
add_executable(fcore_measure fcore_measure.cpp)
target_link_libraries(fcore_measure PRIVATE felitronics::analysis felitronics::oversampling)
target_compile_features(fcore_measure PRIVATE cxx_std_20)

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

# eqviz — renders PNG "quality pictures" (heatmaps + curve sheets) of every eq filter type against
# independent analog references. Offline; catches parameter-corner bugs eyes-first (see tools/eqviz).
add_executable(eqviz eqviz/eqviz.cpp)
target_link_libraries(eqviz PRIVATE felitronics::eq)
target_compile_features(eqviz PRIVATE cxx_std_20)

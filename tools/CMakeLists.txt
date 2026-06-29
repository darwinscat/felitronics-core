# felitronics-core dev tools (built only when this is the top-level project). Not shipped; used to
# validate the DSP against reference implementations (ffmpeg) on real audio.
add_executable(fcore_measure fcore_measure.cpp)
target_link_libraries(fcore_measure PRIVATE felitronics::analysis felitronics::oversampling)
target_compile_features(fcore_measure PRIVATE cxx_std_20)

add_executable(fcore_convolve fcore_convolve.cpp)
target_link_libraries(fcore_convolve PRIVATE felitronics::convolution)
target_compile_features(fcore_convolve PRIVATE cxx_std_20)

# felitronics::io — OFFLINE file I/O: minimal self-contained WAV read/write (16/24/32-bit PCM +
# 32/64-bit float read; 16/24-bit PCM + 32-bit float write; WAVE_FORMAT_EXTENSIBLE; loud rejects).
# Moved verbatim from OrbitCapture's oc/wav.hpp so every product shares one implementation.
# Header-only INTERFACE; zero deps (not even felitronics::core). MESSAGE-THREAD-ONLY (allocates, fopen).
add_library(felitronics_io INTERFACE)
add_library(felitronics::io ALIAS felitronics_io)
target_include_directories(felitronics_io INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(felitronics_io INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_io_tests tests/WavTests.cpp)
    target_link_libraries(felitronics_io_tests PRIVATE felitronics::io)
    target_include_directories(felitronics_io_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_io_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_io_tests COMMAND felitronics_io_tests)

    # The writer's PCM grid is felitronics::dither's: every code survives write∘read, and Dither's output is
    # written as the codes Dither chose. Links dither to TEST that agreement — the io library itself stays
    # zero-dep and restates the rule rather than including it.
    add_executable(felitronics_io_grid_tests tests/WavGridTests.cpp)
    target_link_libraries(felitronics_io_grid_tests PRIVATE felitronics::io felitronics::dither)
    target_include_directories(felitronics_io_grid_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_io_grid_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_io_grid_tests COMMAND felitronics_io_grid_tests)

    # RIFF's word alignment: the pad byte after an odd chunk, written (counted in RIFF, not in the chunk) and
    # skipped on read, with the old unpadded layout and a pad-blind reader planted through the same walk.
    add_executable(felitronics_io_riff_tests tests/WavRiffTests.cpp)
    target_link_libraries(felitronics_io_riff_tests PRIVATE felitronics::io)
    target_include_directories(felitronics_io_riff_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_io_riff_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_io_riff_tests COMMAND felitronics_io_riff_tests)
endif()

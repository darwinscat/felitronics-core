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
endif()

# felitronics::core — Sample alias + size config (SSOT), Math, Smoother, software denormal flush
# (+ optional desktop ScopedFlushToZero). Header-only, zero deps — the shared base.
add_library(felitronics_core INTERFACE)
add_library(felitronics::core ALIAS felitronics_core)
target_include_directories(felitronics_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(felitronics_core INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_core_tests tests/CoreTests.cpp)
    target_link_libraries(felitronics_core_tests PRIVATE felitronics::core)
    target_include_directories(felitronics_core_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_core_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_core_tests COMMAND felitronics_core_tests)

    add_executable(felitronics_core_fft_tests tests/FftTests.cpp)
    target_link_libraries(felitronics_core_fft_tests PRIVATE felitronics::core)
    target_include_directories(felitronics_core_fft_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_core_fft_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_core_fft_tests COMMAND felitronics_core_fft_tests)

    add_executable(felitronics_core_streamresampler_tests tests/StreamResamplerTests.cpp)   # streaming Catmull-Rom rate-match (from OrbitCab)
    target_link_libraries(felitronics_core_streamresampler_tests PRIVATE felitronics::core)
    target_include_directories(felitronics_core_streamresampler_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_core_streamresampler_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_core_streamresampler_tests COMMAND felitronics_core_streamresampler_tests)

    add_executable(felitronics_core_rtstreams_tests tests/RtStreamsTests.cpp)   # RT swap-discipline types (needs threads)
    target_link_libraries(felitronics_core_rtstreams_tests PRIVATE felitronics::core)
    target_include_directories(felitronics_core_rtstreams_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_core_rtstreams_tests PRIVATE cxx_std_20)
    find_package(Threads REQUIRED)
    target_link_libraries(felitronics_core_rtstreams_tests PRIVATE Threads::Threads)
    add_test(NAME felitronics_core_rtstreams_tests COMMAND felitronics_core_rtstreams_tests)
endif()

# felitronics::fftpffft — OPTIONAL, COMPILED SIMD FFT backend (pffft). The FIRST non-header-only module in
# the core: it vendors the 2-file original Julien-Pommier pffft (pffft.c/.h; BSD/FFTPACK — see
# THIRD_PARTY_NOTICES.md) and exposes felitronics::fftpffft::PffftRealFft, a RealFftBackend for the audio
# convolvers. Default OFF (option in the root CMakeLists) so the header-only "link only what you use" and
# "no external deps" contract is unchanged for every consumer that does not opt in.
if(NOT FELITRONICS_WITH_PFFFT)
    return()
endif()

enable_language(C)   # the root project is CXX-only; the vendored pffft.c needs C

# Local-dev guard: under sanitizers the vendored pffft.c and the C++ TUs must share ONE sanitizer runtime. CI
# sets both compilers to clang; a local build with mismatched C/C++ compiler families would link two ASan
# runtimes and fail confusingly. Warn (ci.yml's sanitizer job already passes -DCMAKE_C_COMPILER=clang).
if(FELITRONICS_ENABLE_SANITIZERS AND NOT CMAKE_C_COMPILER_ID STREQUAL CMAKE_CXX_COMPILER_ID)
    message(WARNING "fftpffft + sanitizers: C compiler is ${CMAKE_C_COMPILER_ID}, C++ is ${CMAKE_CXX_COMPILER_ID}. "
                    "Set -DCMAKE_C_COMPILER to the ${CMAKE_CXX_COMPILER_ID} C driver so pffft.c links the same sanitizer runtime.")
endif()

add_library(felitronics_fftpffft STATIC
    src/PffftRealFft.cpp
    pffft/pffft.c)
add_library(felitronics::fftpffft ALIAS felitronics_fftpffft)

target_include_directories(felitronics_fftpffft
    PUBLIC  include
    PRIVATE pffft)                    # vendored <pffft.h> is private to the .cpp — never on a consumer's path
target_link_libraries(felitronics_fftpffft PUBLIC felitronics::core)
target_compile_features(felitronics_fftpffft PUBLIC cxx_std_20)
# Consumers are shared-library plugins (VST3/AU/.so); a non-PIC archive fails to link into a .so on Linux.
# Hidden visibility keeps the vendored pffft's unprefixed C symbols (pffft_*, plus FFTPACK-internal names such
# as cfftf1_ps / validate_pffft_simd) OUT of a consumer plugin's dynamic export table — so two plugins that
# each embed pffft cannot interpose on one another at runtime. (A consumer that ALSO statically links its own
# pffft still gets a loud duplicate-symbol link error — by design; a build should contain a single pffft.)
set_target_properties(felitronics_fftpffft PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_VISIBILITY_PRESET       hidden
    CXX_VISIBILITY_PRESET     hidden
    VISIBILITY_INLINES_HIDDEN ON)
if(UNIX AND NOT APPLE)
    target_link_libraries(felitronics_fftpffft PRIVATE m)   # pffft.c uses libm on Linux
endif()

# The vendored C is warning-clean at -Wall -Wextra, but the DOWNSTREAM-grade set (-Wconversion / -Wcast-align
# / -Wfloat-equal) legitimately fires inside FFT SIMD code. Isolate it: a global -Werror in a parent build
# must never turn pffft.c's inherent casts into errors. Our own adapter .cpp stays under the normal flags.
set_source_files_properties(pffft/pffft.c PROPERTIES COMPILE_OPTIONS
    "$<$<C_COMPILER_ID:GNU,Clang,AppleClang>:-w>$<$<C_COMPILER_ID:MSVC>:/w>")

if(FELITRONICS_BUILD_TESTS)
    add_executable(felitronics_pffft_tests tests/PffftNullTests.cpp)
    target_link_libraries(felitronics_pffft_tests PRIVATE
        felitronics::fftpffft felitronics::convolution felitronics::lineareq)
    target_include_directories(felitronics_pffft_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_pffft_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_pffft_tests COMMAND felitronics_pffft_tests)
endif()

# felitronics::nam — OPTIONAL, COMPILED NeuralAmpModelerCore backend for the process-only neural
# seam. Default OFF: consumers that do not opt in remain header-only and fetch no heavy deps.
if(NOT FELITRONICS_WITH_NAM)
    return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.24)
    message(FATAL_ERROR "FELITRONICS_WITH_NAM requires CMake >= 3.24 ($<LINK_LIBRARY:WHOLE_ARCHIVE> for NAM's architecture self-registration)")
endif()

include(FetchContent)

# NeuralAmpModelerCore ships no usable library target. Populate the pinned sources without entering
# its own project, then compile the same source set used by OrbitCab. Only Eigen's submodule is needed.
set(FELITRONICS_NAM_CORE_TAG "b5a68c3ebed5035a91d9207219346c81e8e3ce8e" CACHE STRING
    "NeuralAmpModelerCore commit to build against")
set(FELITRONICS_NAM_CORE_SRC "" CACHE PATH
    "Local NeuralAmpModelerCore checkout (when empty, fetch the pin)")
if(FELITRONICS_NAM_CORE_SRC AND NOT EXISTS "${FELITRONICS_NAM_CORE_SRC}")
    message(FATAL_ERROR "FELITRONICS_NAM_CORE_SRC does not exist: ${FELITRONICS_NAM_CORE_SRC}")
endif()
if(FELITRONICS_NAM_CORE_SRC)
    message(STATUS "NeuralAmpModelerCore: local checkout — ${FELITRONICS_NAM_CORE_SRC}")
    FetchContent_Declare(nam_core_src
        SOURCE_DIR    "${FELITRONICS_NAM_CORE_SRC}"
        SOURCE_SUBDIR Dependencies/nlohmann)
else()
    message(STATUS "NeuralAmpModelerCore: fetching ${FELITRONICS_NAM_CORE_TAG} from GitHub")
    FetchContent_Declare(nam_core_src
        GIT_REPOSITORY https://github.com/sdatkinson/NeuralAmpModelerCore.git
        GIT_TAG        ${FELITRONICS_NAM_CORE_TAG}
        GIT_SUBMODULES "Dependencies/eigen"
        GIT_SHALLOW    FALSE
        SOURCE_SUBDIR  Dependencies/nlohmann)
endif()
FetchContent_MakeAvailable(nam_core_src)

# namz is a header-only, lossless .nam codec. Pin the v1.1.1 release by immutable commit. SOURCE_SUBDIR
# deliberately points at a directory without CMakeLists.txt: populate the header without building
# namz's project or fetching nlohmann twice.
set(FELITRONICS_NAMZ_TAG "9be9ed8448ab2b72ebfd2605808a4b1f0b24c75a" CACHE STRING
    "namz v1.1.1 commit to build against")
set(FELITRONICS_NAMZ_SRC "" CACHE PATH
    "Local namz checkout (when empty, fetch the pin)")
if(FELITRONICS_NAMZ_SRC AND NOT EXISTS "${FELITRONICS_NAMZ_SRC}")
    message(FATAL_ERROR "FELITRONICS_NAMZ_SRC does not exist: ${FELITRONICS_NAMZ_SRC}")
endif()
if(FELITRONICS_NAMZ_SRC)
    message(STATUS "namz: local checkout — ${FELITRONICS_NAMZ_SRC}")
    FetchContent_Declare(namz_src
        SOURCE_DIR    "${FELITRONICS_NAMZ_SRC}"
        SOURCE_SUBDIR include)
else()
    message(STATUS "namz: fetching ${FELITRONICS_NAMZ_TAG} from GitHub")
    FetchContent_Declare(namz_src
        GIT_REPOSITORY https://github.com/darwinscat/namz.git
        GIT_TAG        ${FELITRONICS_NAMZ_TAG}
        GIT_SHALLOW    TRUE
        SOURCE_SUBDIR  include)
endif()
FetchContent_MakeAvailable(namz_src)

add_library(nam_core STATIC
    ${nam_core_src_SOURCE_DIR}/NAM/activations.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/container.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/conv1d.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/convnet.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/dsp.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/get_dsp.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/linear.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/lstm.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/ring_buffer.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/util.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/wavenet/a2_fast.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/wavenet/model.cpp
    ${nam_core_src_SOURCE_DIR}/NAM/wavenet/slimmable.cpp)
target_compile_features(nam_core PUBLIC cxx_std_20)
target_compile_definitions(nam_core PUBLIC NAM_SAMPLE_FLOAT)
set_target_properties(nam_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(nam_core SYSTEM PUBLIC
    ${nam_core_src_SOURCE_DIR}
    ${nam_core_src_SOURCE_DIR}/NAM
    ${nam_core_src_SOURCE_DIR}/NAM/wavenet
    ${nam_core_src_SOURCE_DIR}/Dependencies/eigen
    ${nam_core_src_SOURCE_DIR}/Dependencies/nlohmann
    ${nam_core_src_SOURCE_DIR}/Dependencies
    ${namz_src_SOURCE_DIR}/include)

add_library(felitronics_nam STATIC src/NamStage.cpp)
add_library(felitronics::nam ALIAS felitronics_nam)
target_include_directories(felitronics_nam PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# NamStage.cpp defines NAMZ_IMPLEMENTATION exactly once. Consumers deliberately receive namz's
# declarations, include <namz.h> without that define, and link the codec symbols from this archive.
# nlohmann rides along as a USAGE REQUIREMENT of those exported headers — namz's rig headers
# (namz_rig*.h) publicly include <nlohmann/json.hpp>, so a consumer TU that includes them must be
# able to resolve it without fetching its own copy. This exports ONLY namz + nlohmann; the NAM/Eigen
# include dirs stay PRIVATE (see the nam_core link below).
target_include_directories(felitronics_nam SYSTEM PUBLIC
    ${namz_src_SOURCE_DIR}/include
    ${nam_core_src_SOURCE_DIR}/Dependencies/nlohmann
    ${nam_core_src_SOURCE_DIR}/Dependencies)
target_compile_features(felitronics_nam PUBLIC cxx_std_20)
target_compile_definitions(felitronics_nam PUBLIC FELITRONICS_WITH_NAM)
set_target_properties(felitronics_nam PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(felitronics_nam
    PRIVATE
    # NAM registers model architectures through static initializers. Encode WHOLE_ARCHIVE here once;
    # plain static linking silently drops those translation units and valid models then fail to load.
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,nam_core>"
    PUBLIC
    felitronics::neural
    felitronics::core)

if(FELITRONICS_BUILD_TESTS)
    # The dedicated RT test controls Eigen's runtime malloc gate from its TU. NAM's Eigen call sites
    # must see the same opt-in macro for the guard to cover allocations originating inside nam_core.
    target_compile_definitions(nam_core PRIVATE EIGEN_RUNTIME_NO_MALLOC)

    add_executable(felitronics_nam_tests tests/NamStageTests.cpp)
    target_link_libraries(felitronics_nam_tests PRIVATE felitronics::nam)
    target_include_directories(felitronics_nam_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
    target_compile_features(felitronics_nam_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_nam_tests COMMAND felitronics_nam_tests)

    add_executable(felitronics_nam_rtalloc_tests tests/NamStageRtAllocTests.cpp)
    target_link_libraries(felitronics_nam_rtalloc_tests PRIVATE felitronics::nam)
    target_include_directories(felitronics_nam_rtalloc_tests PRIVATE
        ${felitronics_core_SOURCE_DIR}/test_support)
    target_include_directories(felitronics_nam_rtalloc_tests SYSTEM PRIVATE
        ${nam_core_src_SOURCE_DIR}/Dependencies/eigen)
    target_compile_definitions(felitronics_nam_rtalloc_tests PRIVATE EIGEN_RUNTIME_NO_MALLOC)
    target_compile_features(felitronics_nam_rtalloc_tests PRIVATE cxx_std_20)
    add_test(NAME felitronics_nam_rtalloc_tests COMMAND felitronics_nam_rtalloc_tests)
endif()

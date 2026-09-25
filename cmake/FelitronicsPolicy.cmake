# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#
# felitronics build POLICY — the compile/link flags every tree that compiles felitronics code states for
# itself: FP contraction (law 10), the MSVC stack and UTF-8, the emscripten physics, the sanitizer option,
# and the strict warning set of the header-hygiene gate.
#
# include()d, not add_subdirectory()'d: add_compile_options() applies to the directory that includes this
# file and everything below it, and to nothing outside. So core includes it from its top-level
# CMakeLists.txt, and a satellite repository (felitronics-guitar-core) includes the SAME file from its
# own — a consumer that FetchContent's core does not inherit these flags for its own targets, and a
# satellite that restated them by hand would be a second copy of law 10.
#   include(${felitronics_core_SOURCE_DIR}/cmake/FelitronicsPolicy.cmake)
if(MSVC)
    add_compile_options(/utf-8)         # sources use UTF-8 in comments (—, φ, ≈, ·) — MSVC needs this flag
    add_link_options(/STACK:8388608)    # 8 MB stack (the Unix default). Tests put a ~200 KB EqEngine on the
                                        # stack; MSVC's 1 MB default overflows → SegFault. Match Unix.
else()
    # FP CONTRACTION, stated instead of inherited (DSP-ARCHITECTURE.md §2 law 10). `a*b + c` may become a
    # single FMA with ONE rounding instead of two, and the toolchains do not agree on when: clang defaults
    # to `on` (fusion within one expression, the C/C++ standard's own rule), GCC defaults to `fast` (fusion
    # ACROSS statements). That difference is not academic — it broke two exactness claims the moment an
    # arm64 Linux row reached CI, on untouched main: a nonlinear stage's "silence in ⇒ output identically 0" read
    # 1.26e-08, and the multi-res fast pane missed its budget by 0.005841 dB. Neither the Mac (clang) nor
    # any x86 row could see it: baseline x86-64 has no FMA instruction to contract with at all.
    #
    # `on`, not `off`. It fixes both (verified on gcc 13.3 and 14.2, aarch64) and costs the shipping tier
    # NOTHING, because it is already clang's default. `off` also fixes them but disables FMA everywhere,
    # measured at +12% on the Svf sample loop (16.25 -> 18.24 ms, arm64 macOS) for no demonstrated benefit.
    # NB on GCC <= 13 `on` is implemented as `off` (no fmadd emitted at all); GCC 14 emits the standard
    # behaviour — 1 fmadd for the single-expression form, 0 for the cross-statement one. Both pass.
    # MSVC needs no equivalent: its default /fp:precise does not contract.
    add_compile_options(-ffp-contract=on)
endif()

# The `wasm-audio` tier's PHYSICS. Its POLICY (-fno-exceptions -fno-rtti) lives in the `wasm-audio` CMake
# preset; these three are needed by ANY emscripten build of this tree, so they are not left to a preset a
# developer might forget. All three were measured on emsdk 6.0.9, not guessed:
#   STACK_SIZE           emscripten defaults to 64 KB. The suite puts a ~200 KB EqEngine on the stack (the
#                        same fact the MSVC /STACK line above exists for), and in wasm a blown stack does NOT
#                        reliably trap: it first produced WRONG ANSWERS in the analyser tap tests (a reported
#                        hop of 600 against a requested 512) and only then "memory access out of bounds".
#                        8 MB = the Unix default the MSVC line already matches.
#   ALLOW_MEMORY_GROWTH  the default heap is a FIXED 16 MB; the NUPC stress suite and two measurement suites
#                        legitimately need more and died with Aborted(OOM).
#   ASSERTIONS           without it a failure prints a bare "Aborted()" and nothing else. A gate whose red
#                        light says nothing is a bad gate, and these binaries are never shipped. It is a
#                        cache variable because a directory-level link option lands AFTER
#                        CMAKE_EXE_LINKER_FLAGS, so a hard-coded =1 here would silently DOWNGRADE the
#                        checked configuration's -sASSERTIONS=2 (emscripten keeps the last spelling —
#                        verified in the generated link.txt). The checked CI step passes 2 through this.
set(FELITRONICS_WASM_ASSERTIONS "1" CACHE STRING
    "-sASSERTIONS level for emscripten binaries (2 in the checked configuration)")
if(EMSCRIPTEN)
    add_link_options(-sSTACK_SIZE=8388608 -sALLOW_MEMORY_GROWTH=1 -sASSERTIONS=${FELITRONICS_WASM_ASSERTIONS})
    # ...and the one flag that makes "adding a thread fails the build" TRUE rather than aspirational.
    # Without it emscripten links pthread STUBS whose pthread_create returns ENOTSUP, so std::thread
    # compiles, links, and only aborts at runtime — measured. --wrap turns the reference into an undefined
    # __wrap_pthread_create at link time: `wasm-ld: error: undefined symbol: __wrap_pthread_create`. It
    # works because libc++'s thread constructor is header-inline, so the reference lands in OUR object.
    # Measured zero false positives across all 78 links in this tree.
    # It does NOT catch std::mutex — those live inside libc++-noexcept.a and the stub silently succeeds —
    # so the CI job also lints module headers for <thread>/<mutex>/<condition_variable>/<future>.
    # Deliberately on for EVERY emscripten build of this tree, not only the preset: DSP-ARCHITECTURE.md 2
    # declares `wasm-audio` single-threaded, so a threaded wasm build of this core is not a thing we mean to
    # allow silently. If that ever changes it is one line, changed on purpose.
    add_link_options(LINKER:--wrap=pthread_create)
endif()

# ASan + UBSan for dev + a dedicated CI job. UBSan makes platform-divergent UB LOUD on any host: a
# convolution divide-by-zero on an unprepared engine was SILENT on Apple-Silicon (ARM SDIV /0 → 0) yet
# crashed on x86-64 (IDIV #DE) and would HardFault on bare-metal. `-fno-sanitize-recover=all` turns every
# UB into a hard failure (== a trap), so a future port to "who knows what" (Cortex-M …) stays honest.
# Flags are global so every module self-test is instrumented. Clang/GCC only (MSVC has no UBSan).
option(FELITRONICS_ENABLE_SANITIZERS "Build the self-tests with ASan + UBSan (hard-fail; Clang/GCC)" OFF)
if(FELITRONICS_ENABLE_SANITIZERS)
    if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU"))
        message(FATAL_ERROR "FELITRONICS_ENABLE_SANITIZERS requires Clang or GCC (MSVC has no UBSan).")
    endif()
    message(STATUS "felitronics-core: ASan + UBSan enabled (hard-fail on any UB)")
    add_compile_options(-g -fno-omit-frame-pointer
                        -fsanitize=address,undefined -fno-sanitize-recover=all)  # NOT -fsanitize=integer: the dither RNG wraps unsigned on purpose
    add_link_options(-fsanitize=address,undefined)
    add_compile_definitions(_GLIBCXX_ASSERTIONS)                                  # libstdc++ container bounds checks (inert on libc++; ASan covers OOB there)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-ftrivial-auto-var-init=pattern)                     # cheap uninitialised-stack tripwire
    endif()
endif()

# The header-hygiene gate's warning set (Clang/GCC). One list, used by core's hygiene TU and by any
# satellite's, so "clean under the gate" means the same thing in every tree.
set(FELITRONICS_HYGIENE_WARNINGS
    -Werror -Wall -Wextra -Wpedantic
    -Wconversion -Wsign-conversion -Wfloat-equal -Wshadow -Wcast-align
    -Wzero-as-null-pointer-constant -Wunreachable-code -Wswitch-enum -Wextra-semi
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wshorten-64-to-32;-Wimplicit-float-conversion;-Wconditional-uninitialized;-Wconstant-conversion;-Wbool-conversion>)

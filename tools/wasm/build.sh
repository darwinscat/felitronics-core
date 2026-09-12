#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#
# Builds the wasm artifacts of this repo's two C ABIs:
#
#   fcprobe.*   — the P0 measurement spike (tools/wasm/fc_probe.cpp)
#   fcmaster.*  — the mastering ABI (tools/fc_master_abi.h + tools/wasm/fc_master.cpp), which is what a
#                 browser worker links to render and to run the loudness search. Same flags, same numeric
#                 contract, same no-threads and byte-identity checks — see the fc_master section below.
#
# The probe part first. Three artifacts from one source:
#
#   fcprobe.web.js/.wasm    -sENVIRONMENT=web,worker — EXACTLY what P0 specifies. This is the artifact whose
#                           size is reported and which the page loads; the no-threads claim is made about it.
#   fcprobe.node.js/.wasm   -sENVIRONMENT=node — the same wasm with node glue, so the parity harness runs
#                           without a browser (and so CI can). The script verifies the two .wasm files are
#                           BYTE-IDENTICAL, which is what lets a parity result proven on one transfer to the
#                           other: -sENVIRONMENT only shapes JS glue.
#   fcprobe.debug.js/.wasm  -O1 -g -sASSERTIONS=2 -sSAFE_HEAP=1 -sSTACK_OVERFLOW_CHECK=2 — the libsoxr lesson.
#                           libsoxr compiled and linked cheaply under emscripten and crashed at RUNTIME; a
#                           release build that "works" proves less than a checked build that agrees with it.
#
# Requires emsdk on PATH:  source ~/emsdk/emsdk_env.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="${1:-$HERE/build}"
mkdir -p "$OUT"

command -v em++ >/dev/null || { echo "em++ not on PATH — source \$EMSDK/emsdk_env.sh first"; exit 1; }
echo "emcc: $(emcc --version | head -1)"

INC=(-I"$ROOT/tools"
     -I"$ROOT/modules/core/include"
     -I"$ROOT/modules/analysis/include"
     -I"$ROOT/modules/oversampling/include")

# The mastering ABI pulls in the whole chain. This list is `felitronics::mastering`'s own link list in
# modules/mastering/CMakeLists.txt, spelled as include paths — plus oversampling, which analysis needs.
# There is nothing to compile but the two headers' worth of templates: every one of these modules is
# header-only (INTERFACE libraries), which is why one em++ invocation is the whole build.
MASTER_INC=(-I"$ROOT/tools"
            -I"$ROOT/modules/mastering/include"
            -I"$ROOT/modules/core/include"
            -I"$ROOT/modules/eq/include"
            -I"$ROOT/modules/dynamics/include"
            -I"$ROOT/modules/saturation/include"
            -I"$ROOT/modules/limiter/include"
            -I"$ROOT/modules/stereo/include"
            -I"$ROOT/modules/dither/include"
            -I"$ROOT/modules/analysis/include"
            -I"$ROOT/modules/oversampling/include")

# The numeric contract. -ffp-contract=off is stated explicitly on BOTH sides rather than relied on: baseline
# wasm has no scalar FMA so emscripten cannot contract anyway, but saying so keeps the two build files
# symmetric and survives a future toolchain that grows the ability.
# NOT -mrelaxed-simd, ever: f64x2.relaxed_madd is implementation-defined (fused on hosts with FMA, unfused
# elsewhere), which breaks determinism between MACHINES, not merely between tiers.
NUMERIC=(-ffp-contract=off -fno-fast-math)

COMMON=(-std=c++20 -fno-exceptions -fno-rtti "${NUMERIC[@]}" "${INC[@]}"
        --no-entry
        -sMODULARIZE=1
        -sEXPORT_NAME=createFcProbe
        -sALLOW_MEMORY_GROWTH=1
        -sFILESYSTEM=0
        -sMALLOC=emmalloc
        "-sEXPORTED_FUNCTIONS=['_fc_probe_run','_fc_probe_lufs','_fc_probe_dbtp','_fc_probe_tp_linear','_fc_probe_sample_peak','_fc_probe_block_count','_fc_probe_block_energies','_fc_probe_dropped','_fc_probe_os_factor','_fc_probe_os_taps','_fc_probe_chunk','_fc_probe_sizeof_longdouble','_malloc','_free']"
        "-sEXPORTED_RUNTIME_METHODS=['HEAPF32','HEAPF64']")
        # _malloc/_free and the HEAP views are OPT-IN in emscripten 6.x — without these two lines
        # Module._malloc and Module.HEAPF32 are simply `undefined` and the page dies on first use.

SRC="$HERE/fc_probe.cpp"

echo "--- web (the P0 artifact)"
em++ "${COMMON[@]}" -O3 -sENVIRONMENT=web,worker "$SRC" -o "$OUT/fcprobe.web.js"

echo "--- node (same wasm, node glue — for the parity harness)"
em++ "${COMMON[@]}" -O3 -sENVIRONMENT=node "$SRC" -o "$OUT/fcprobe.node.js"

echo "--- debug (SAFE_HEAP + assertions + stack checks)"
em++ "${COMMON[@]}" -O1 -g -sASSERTIONS=2 -sSAFE_HEAP=1 -sSTACK_OVERFLOW_CHECK=2 \
     -sENVIRONMENT=node "$SRC" -o "$OUT/fcprobe.debug.js"

echo
echo "=== the two release .wasm must be byte-identical (ENVIRONMENT shapes glue, not code)"
a=$(shasum -a 256 "$OUT/fcprobe.web.wasm"  | cut -d' ' -f1)
b=$(shasum -a 256 "$OUT/fcprobe.node.wasm" | cut -d' ' -f1)
echo "  web  $a"
echo "  node $b"
[ "$a" = "$b" ] && echo "  IDENTICAL" || { echo "  *** DIFFER — a parity result on one does not transfer to the other"; exit 1; }

echo
echo "=== no threads, proven from the artifact rather than from the page loading"
# Parses the wasm memory section directly, so this checks something on every machine — the earlier
# wasm-objdump step silently checked NOTHING wherever that tool was not installed, which is most machines.
node "$HERE/check-no-threads.mjs" "$OUT/fcprobe.web.wasm" "$OUT/fcprobe.web.js"

echo
echo "=== size (acceptance criterion 2)"
printf "  %-22s %10s %10s %10s\n" file raw gzip brotli
for f in fcprobe.web.wasm fcprobe.web.js; do
    raw=$(wc -c < "$OUT/$f")
    gz=$(gzip -9 -c "$OUT/$f" | wc -c)
    br=$(brotli -q 11 -c "$OUT/$f" 2>/dev/null | wc -c || echo "n/a")
    printf "  %-22s %10s %10s %10s\n" "$f" "$raw" "$gz" "$br"
done

#==================================================================================================
# fc_master — the MASTERING ABI (tools/fc_master_abi.h, implemented by tools/wasm/fc_master.cpp)
#
# Built exactly like fc_probe above and for the same reasons; only three things differ, and each is a
# consequence of what this ABI is rather than a preference:
#
#  1. THE EXPORT LIST IS GENERATED FROM THE SOURCE, not typed here. `-sEXPORTED_FUNCTIONS` is a
#     whitelist: a name missing from it is dead-stripped, and the page then finds `Module._fc_…`
#     undefined at the moment it needs it — at runtime, in a worker, with no build-time diagnostic
#     anywhere. EMSCRIPTEN_KEEPALIVE alone does NOT save it once EXPORTED_FUNCTIONS is given. So the
#     list is read out of the FC_EXPORT lines of fc_master.cpp: an entry point added to the ABI is
#     exported by the fact of existing, and cannot be forgotten here.
#  2. -sSTACK_SIZE=8388608. The repository's CMakeLists gives this to every emscripten build of this
#     tree and says why (a blown wasm stack does not reliably trap — it produced WRONG ANSWERS before
#     it produced an out-of-bounds). `MasteringChain` heap-allocates, but the solver's per-pass work
#     and this file's own planar pointer tables are automatic, and 64 KB — emscripten's default — is
#     not a number anybody chose for them.
#  3. -msimd128, which is NOT a hopeful flag here. Since v0.31.0 the three polyphase loops of this tree
#     are one function, `core::firDot`, with four hand-written kernels; the wasm one is behind
#     `__wasm_simd128__` and is simply absent without this flag. Measured on a 5:21 stereo programme,
#     one render pass: 8748 ms at v0.30.0 -> 4339 ms on v0.31.0's repacking alone -> **2802 ms** with the
#     kernel, i.e. 3.1x, and a full solve 29.3 s -> 10.6 s.
#     IT IS NUMERICALLY FREE, and that is measured rather than taken on trust: the same programme
#     through the scalar and the SIMD module agrees in ALL 30 877 716 delivered samples BIT FOR BIT, and
#     so do the gain, the ceiling, the integrated loudness and the true peak to 17 digits. That is the
#     kernels' own design — four partial accumulators, (s0+s1)+(s2+s3), every operation rounded
#     separately, identical in all four — and the core gates it with
#     `felitronics_core_polyphasefir_tests`. Which is why it goes on THIS target and not on fc_probe:
#     the probe's acceptance is its own, and a flag is not added to a contract it was not measured on.
#  4. -sEXPORT_NAME=createFcMaster, so a page that loads BOTH modules gets two factories rather than
#     one name overwriting the other. Everything else — the numeric contract, emmalloc, no filesystem,
#     growable memory, the HEAPF32/HEAPF64 opt-in — is the probe's line for the probe's reasons.
#
# NOT built: a debug variant. fc_probe has one because libsoxr taught that a release build that "works"
# proves less than a checked one that agrees with it; the checked configuration of THIS ABI already
# exists and is stronger — `felitronics_master_abi_tests`, `felitronics_master_abi_reentry_tests` and
# `fcore_master selftest` run the same fc_master.cpp under ctest with ASan and UBSan.
#==================================================================================================
MSRC="$HERE/fc_master.cpp"

# The whitelist, read from the source. `-o` prints one name per match, so an entry point that shares a
# line with another is still found. _malloc/_free are the page's own, and are opt-in in emscripten 6.x.
MEXPORTS=$(grep -oE 'FC_EXPORT[[:space:]]+(void|std::uint32_t|uint32_t|fc_status)[[:space:]]+fc_[a-z_]+' "$MSRC" \
           | awk '{print "_" $NF}' | sort -u | paste -sd, -)
MEXPORTS="$MEXPORTS,_malloc,_free"
echo "--- fc_master exports: $(printf '%s\n' "$MEXPORTS" | tr ',' '\n' | wc -l | tr -d ' ') symbols"
# A generated list that silently came back empty would build a module with nothing in it, so it is
# checked rather than trusted. 25 entry points at ABI v1; the guard is a floor, not the count, so
# appending one is not a build break.
[ "$(printf '%s\n' "$MEXPORTS" | tr ',' '\n' | wc -l | tr -d ' ')" -ge 27 ] \
    || { echo "*** the export list did not come out of $MSRC — refusing to link a module with no ABI"; exit 1; }

MCOMMON=(-std=c++20 -fno-exceptions -fno-rtti "${NUMERIC[@]}" "${MASTER_INC[@]}"
         -msimd128
         --no-entry
         -sMODULARIZE=1
         -sEXPORT_NAME=createFcMaster
         -sALLOW_MEMORY_GROWTH=1
         -sSTACK_SIZE=8388608
         -sFILESYSTEM=0
         -sMALLOC=emmalloc
         "-sEXPORTED_FUNCTIONS=[$MEXPORTS]"
         "-sEXPORTED_RUNTIME_METHODS=['HEAPF32','HEAPF64']")

# -sEXPORT_ES6 on the WEB variant only. The page loads this from a module worker, where an ES module
# with a default export is the only shape an `import` can take; node's harness here keeps the plain
# MODULARIZE factory so `createRequire` still works. It shapes GLUE and nothing else — the byte-identity
# check below is what says so rather than the flag's documentation.
echo "--- fc_master web (ES module, for a module worker)"
em++ "${MCOMMON[@]}" -O3 -sENVIRONMENT=web,worker -sEXPORT_ES6=1 "$MSRC" -o "$OUT/fcmaster.web.mjs"

echo "--- fc_master node (same wasm, node glue — for the parity harness)"
em++ "${MCOMMON[@]}" -O3 -sENVIRONMENT=node "$MSRC" -o "$OUT/fcmaster.node.js"

echo
echo "=== the two fc_master .wasm must be byte-identical, same rule as above"
a=$(shasum -a 256 "$OUT/fcmaster.web.wasm"  | cut -d' ' -f1)
b=$(shasum -a 256 "$OUT/fcmaster.node.wasm" | cut -d' ' -f1)
echo "  web  $a"
echo "  node $b"
[ "$a" = "$b" ] && echo "  IDENTICAL" || { echo "  *** DIFFER"; exit 1; }

echo
echo "=== no threads (fc_master)"
node "$HERE/check-no-threads.mjs" "$OUT/fcmaster.web.wasm" "$OUT/fcmaster.web.mjs"

echo
echo "=== size"
printf "  %-22s %10s %10s %10s\n" file raw gzip brotli
for f in fcmaster.web.wasm fcmaster.web.mjs; do
    raw=$(wc -c < "$OUT/$f")
    gz=$(gzip -9 -c "$OUT/$f" | wc -c)
    br=$(brotli -q 11 -c "$OUT/$f" 2>/dev/null | wc -c || echo "n/a")
    printf "  %-22s %10s %10s %10s\n" "$f" "$raw" "$gz" "$br"
done

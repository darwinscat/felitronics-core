#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#
# Builds the P0 wasm spike. Three artifacts from one source:
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

#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Validate felitronics-core analysis (LUFS, true-peak) against ffmpeg's ebur128 reference, on REAL songs
# and on SYNTHETIC adversarial cases chosen to expose classic implementation bugs. ffmpeg and our tool
# are always measured on the SAME decoded signal (mono↔stereo up-mix applies a -3 dB pan law — comparing
# different signals is itself a trap).
#
# Usage: tools/validate_ffmpeg.sh [build-dir]   (run from the repo root)
set -u
BUILD="${1:-build}"
M="$BUILD/tools/fcore_measure"
AUDIO="/Users/oleh/Downloads/AUDIO"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
command -v ffmpeg >/dev/null || { echo "no ffmpeg"; exit 1; }
[ -x "$M" ] || { echo "build $M first"; exit 1; }

refI()  { ffmpeg -nostats -i "$1" -af ebur128 -f null - 2>&1 | grep -oE 'I:[[:space:]]+-?[0-9.]+' | tail -1 | grep -oE '\-?[0-9.]+'; }
refTP() { ffmpeg -nostats -i "$1" -af 'ebur128=peak=true' -f null - 2>&1 | grep -A1 'True peak' | grep -oE 'Peak:[[:space:]]+-?[0-9.]+' | grep -oE '\-?[0-9.]+'; }
d()     { awk -v x="$1" -v y="$2" 'BEGIN{printf "%+.2f", x-y}'; }

printf "%-34s | %-26s | %-26s\n" "case" "LUFS  ref / mine / Δ" "TruePeak  ref / mine / Δ"
printf -- "----------------------------------+----------------------------+----------------------------\n"

check() { # name src sr ch
  local name="$1" src="$2" sr="$3" ch="$4"
  ffmpeg -v error -y -i "$src" -map a:0 -f f32le "$TMP/x.raw" 2>/dev/null
  local rI rT mI mT
  rI=$(refI "$src"); rT=$(refTP "$src")
  mI=$("$M" lufs "$sr" "$ch" "$TMP/x.raw"); mT=$("$M" truepeak "$sr" "$ch" "$TMP/x.raw")
  printf "%-34s | %7s /%7s /%6s | %7s /%7s /%6s\n" "$name" "$rI" "$mI" "$(d "$mI" "$rI")" "$rT" "$mT" "$(d "$mT" "$rT")"
}

echo "### REAL SONGS"
for f in cat-in-space-2020_DEAR42052542 cats-evolution-DEAR42037456 cat-in-space-2025_QT6G32570545; do
  src="$AUDIO/$f.flac"; [ -f "$src" ] || continue
  read sr ch < <(ffprobe -v error -select_streams a:0 -show_entries stream=sample_rate,channels -of csv=p=0 "$src" | tr ',' ' ')
  check "$f" "$src" "$sr" "$ch"
done

echo "### SYNTHETIC adversarial (expose: ISP, gating, coeff-recompute)"
gen() { ffmpeg -v error -y -f lavfi -i "$1" "$2" 2>/dev/null; }
# ISP: 12 kHz sine phased to miss the sample crests → true peak >> sample peak (sample-peak meters under-read)
gen "aevalsrc=0.9*sin(2*PI*12000*t+0.7)|0.9*sin(2*PI*12000*t+0.7):s=48000:d=5" "$TMP/isp.wav";  check "ISP 12kHz phased (48k)" "$TMP/isp.wav" 48000 2
# coeff recompute: 1 kHz at 44.1 k and 96 k (a hardcoded-48k K-weighting would diverge here)
gen "aevalsrc=0.3*sin(2*PI*1000*t)|0.3*sin(2*PI*1000*t):s=44100:d=10" "$TMP/s441.wav"; check "1kHz @44.1k (coeff recompute)" "$TMP/s441.wav" 44100 2
gen "aevalsrc=0.3*sin(2*PI*1000*t)|0.3*sin(2*PI*1000*t):s=96000:d=10" "$TMP/s96.wav";  check "1kHz @96k  (coeff recompute)" "$TMP/s96.wav"  96000 2
# gating: 5 s loud then 10 s near-silent → integrated must reflect the LOUD part (relative gate excludes quiet)
gen "aevalsrc='0.5*sin(2*PI*1000*t)*if(lt(t\,5)\,1\,0.004)'|'0.5*sin(2*PI*1000*t)*if(lt(t\,5)\,1\,0.004)':s=48000:d=15" "$TMP/ql.wav"; check "loud 5s + quiet 10s (gate)" "$TMP/ql.wav" 48000 2
# near-silence: everything below the -70 absolute gate
gen "aevalsrc=0.002*sin(2*PI*1000*t)|0.002*sin(2*PI*1000*t):s=48000:d=5" "$TMP/q.wav"; check "near-silence (-70 abs gate)" "$TMP/q.wav" 48000 2
# full-scale square-ish (broadband ISP + the hardest true-peak): a 0 dBFS sine near fs/4
gen "aevalsrc=0.999*sin(2*PI*11000*t+1.1)|0.999*sin(2*PI*11000*t+1.1):s=48000:d=5" "$TMP/fs.wav"; check "near-FS 11kHz (true-peak)" "$TMP/fs.wav" 48000 2

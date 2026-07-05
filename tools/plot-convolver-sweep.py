#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#
# Render docs/assets/nupc-vs-juce.svg — the two-machine NUPC-vs-juce::dsp::Convolution DAW-buffer sweep — from
# fcore_fftbench CSV. Self-contained SVG (inline presentation attributes, its own light background) so GitHub
# renders it via <img> in markdown on both light and dark themes.
#
# Produce the CSVs (one per machine) with the FINE log-ladder sweep, then render:
#   FCORE_FINE_SWEEP=1 ./build-juce/tools/fcore_fftbench | grep -E '^  [0-9]' > i9.csv   # on the i9 (pin to a P-core)
#   FCORE_FINE_SWEEP=1 ./build-juce/tools/fcore_fftbench | grep -E '^  [0-9]' > m5.csv   # on the M5
#   python3 tools/plot-convolver-sweep.py i9.csv m5.csv docs/assets/nupc-vs-juce.svg
# Fine CSV columns: buffer[*],ms,juce_mean,juce_spread%,nupc_mean,nupc_spread%,nupc_worst (juce=-1 below 16, where it
# is not measured). Each %RT is the MEAN of 3–10 warmed windows (adaptive). Legacy 5-col CSV is also accepted.
import math, re, sys

def load(p):
    B, j, u, m = [], [], [], []
    for line in open(p):
        parts = [x.strip() for x in line.split(',')]
        if len(parts) < 5: continue
        try: b = int(re.sub(r'[^0-9]', '', parts[0]))
        except ValueError: continue
        if len(parts) >= 7:                       # fine: buffer,ms,juce_mean,juce_spread,nupc_mean,nupc_spread,nupc_worst
            jv, uv, mv = float(parts[2]), float(parts[4]), float(parts[6])
        else:                                     # legacy: buffer,ms,juce,nupc,nupc_max
            jv, uv, mv = float(parts[2]), float(parts[3]), float(parts[4])
        B.append(b); j.append(jv if jv >= 0 else None); u.append(uv); m.append(mv)
    return B, j, u, m

i9_csv, m5_csv, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
B, ji, ui, mi = load(i9_csv)
_, jm, um, mm = load(m5_csv)

INK, INK2, INK3 = "#1a2230", "#55606b", "#8a95a0"
GRID, BORDER, PANEL = "#e6ebf0", "#dde3e9", "#ffffff"
US, JUCE, DANGER = "#2a78d6", "#eb6834", "#d03b3b"
MONO = "ui-monospace,'SF Mono',Menlo,Consolas,monospace"
SANS = "system-ui,-apple-system,'Segoe UI',Roboto,sans-serif"
W = 880
esc = lambda s: str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
pow2 = lambda b: (b & (b - 1)) == 0

def cross_idx(j, u):
    n = len(j)
    for i in range(n - 1):
        if j[i] is not None and j[i + 1] is not None and j[i] < u[i] and j[i + 1] < u[i + 1]: return i
    return n - 1 if (j[-1] is not None and j[-1] < u[-1]) else n

def panel(x0, y0, cpu, env, juce, us, usMax):
    L, R, T = 52, 108, 66
    pw, ph = W - L - R, 300
    x1, y1 = x0 + L, y0 + T
    lb0, lb1 = math.log2(B[0]), math.log2(B[-1])
    yLo, yHi = 0.15, 300.0
    def X(i): return x1 + pw * (math.log2(B[i]) - lb0) / (lb1 - lb0)
    def Y(v):
        v = max(v, yLo)
        return y1 + ph * (1 - (math.log10(v) - math.log10(yLo)) / (math.log10(yHi) - math.log10(yLo)))
    s = []
    jvals = [v for v in juce if v is not None]
    j0 = next(i for i in range(len(B)) if juce[i] is not None)   # first measured JUCE point (16)
    usMean = sum(us) / len(us); jMin = min(jvals)
    ci = cross_idx(juce, us); crossNext = B[ci] if ci < len(B) else B[-1]
    s.append(f'<text x="{x0+L}" y="{y0+22}" font-family="{SANS}" font-size="16" font-weight="700" fill="{INK}">{esc(cpu)}</text>')
    s.append(f'<text x="{x0+L}" y="{y0+40}" font-family="{MONO}" font-size="11.5" fill="{INK3}">{esc(env)}</text>')
    cap = (f'NUPC flat mean {usMean:.2f}%  ·  JUCE {jMin:.2f}–{juce[j0]:.0f}%  ·  '
           f'NUPC cheaper up to ~{crossNext} samples  ·  0 ms latency')
    s.append(f'<text x="{x0+L}" y="{y0+56}" font-family="{MONO}" font-size="11.5" fill="{INK2}">{esc(cap)}</text>')
    if 0 < ci < len(B):
        xc = (X(ci - 1) + X(ci)) / 2
        s.append(f'<rect x="{x1:.1f}" y="{y1:.1f}" width="{xc-x1:.1f}" height="{ph}" fill="{US}" opacity="0.055"/>')
        s.append(f'<line x1="{xc:.1f}" y1="{y1:.1f}" x2="{xc:.1f}" y2="{y1+ph:.1f}" stroke="{INK3}" stroke-width="1" stroke-dasharray="4 4" opacity="0.5"/>')
        s.append(f'<text x="{x1+6:.1f}" y="{y1+14:.1f}" font-family="{MONO}" font-size="11" fill="{US}">◄ NUPC cheaper (mean)</text>')
        if xc < x1 + pw - 120:
            s.append(f'<text x="{xc+7:.1f}" y="{y1+14:.1f}" font-family="{MONO}" font-size="11" fill="{JUCE}">JUCE cheaper ►</text>')
    for g in [0.3, 1, 3, 10, 30, 100, 300]:
        yy = Y(g); ceil = g == 100
        dash = ' stroke-dasharray="6 4"' if ceil else ''
        s.append(f'<line x1="{x1:.1f}" y1="{yy:.1f}" x2="{x1+pw:.1f}" y2="{yy:.1f}" stroke="{DANGER if ceil else GRID}" stroke-width="{1.4 if ceil else 1}"{dash}/>')
        lab = f"{g:.1f}%" if g < 1 else f"{g}%"
        s.append(f'<text x="{x1-8:.1f}" y="{yy+4:.1f}" font-family="{MONO}" font-size="11" fill="{INK3}" text-anchor="end">{lab}</text>')
    s.append(f'<text x="{x1+pw-2:.1f}" y="{Y(100)-5:.1f}" font-family="{MONO}" font-size="11" fill="{DANGER}" text-anchor="end">100% real-time</text>')
    for i, b in enumerate(B):
        if not pow2(b): continue
        s.append(f'<text x="{X(i):.1f}" y="{y1+ph+19:.1f}" font-family="{MONO}" font-size="11" fill="{INK3}" text-anchor="middle">{b}</text>')
    s.append(f'<text x="{x1+pw/2:.1f}" y="{y1+ph+38:.1f}" font-family="{MONO}" font-size="11" fill="{INK2}" text-anchor="middle">host buffer size (samples, log scale · &lt;16 = synthetic, not a real host buffer)</text>')
    # worst-buffer band (nupc mean -> worst)
    pts = [f"{X(i):.1f} {Y(us[i]):.1f}" for i in range(len(B))]
    pts += [f"{X(i):.1f} {Y(usMax[i]):.1f}" for i in range(len(B) - 1, -1, -1)]
    s.append(f'<path d="M{" L".join(pts)} Z" fill="{US}" opacity="0.10"/>')
    # lines (JUCE only where measured, i.e. b>=16)
    jpts = [f"{X(i):.1f} {Y(juce[i]):.1f}" for i in range(len(B)) if juce[i] is not None]
    s.append(f'<path d="M{" L".join(jpts)}" fill="none" stroke="{JUCE}" stroke-width="2" stroke-linejoin="round"/>')
    upts = [f"{X(i):.1f} {Y(us[i]):.1f}" for i in range(len(B))]
    s.append(f'<path d="M{" L".join(upts)}" fill="none" stroke="{US}" stroke-width="2" stroke-linejoin="round"/>')
    li = len(B) - 1
    s.append(f'<text x="{X(li)+9:.1f}" y="{Y(us[li])+4:.1f}" font-family="{MONO}" font-size="12" font-weight="600" fill="{US}">NUPC {us[li]:.1f}%</text>')
    s.append(f'<text x="{X(li)+9:.1f}" y="{Y(juce[li])+4:.1f}" font-family="{MONO}" font-size="12" font-weight="600" fill="{JUCE}">JUCE {juce[li]:.1f}%</text>')
    s.append(f'<text x="{X(j0)+4:.1f}" y="{Y(juce[j0])-7:.1f}" font-family="{MONO}" font-size="12" font-weight="600" fill="{JUCE}">{juce[j0]:.0f}%</text>')
    return "\n".join(s)

PH = 412
H = 92 + 2 * PH
out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="{SANS}">',
       f'<rect x="0" y="0" width="{W}" height="{H}" rx="14" fill="{PANEL}" stroke="{BORDER}"/>',
       f'<text x="28" y="34" font-family="{SANS}" font-size="19" font-weight="700" fill="{INK}">Convolver CPU across the full DAW buffer range</text>',
       f'<text x="28" y="55" font-family="{MONO}" font-size="11.5" fill="{INK3}">MatrixConvolverNupc vs juce::dsp::Convolution 8.0.4 · 131072-tap stereo (LRDiag) · 48 kHz · mean of 3–10 warmed windows</text>']
lx = 28
for col, txt, band in [(US, "NUPC mean", False), (US, "NUPC worst buffer (band)", True), (JUCE, "juce::dsp::Convolution (from 16)", False), (DANGER, "100% real-time", False)]:
    dash = ' stroke-dasharray="5 3"' if col == DANGER else ''
    out.append(f'<line x1="{lx}" y1="76" x2="{lx+22}" y2="76" stroke="{col}" stroke-width="{8 if band else 2.6}" opacity="{0.22 if band else 1}"{dash}/>')
    out.append(f'<text x="{lx+28}" y="80" font-family="{SANS}" font-size="12" fill="{INK2}">{esc(txt)}</text>')
    lx += 32 + len(txt) * 7.0 + 24
out.append(panel(0, 92, "Intel i9-13900H", "Debian 13 · gcc 14 · pffft/SSE · P-core pinned", ji, ui, mi))
out.append(panel(0, 92 + PH, "Apple M5 Pro", "macOS · pffft/NEON · P-core", jm, um, mm))
out.append('</svg>')
open(out_path, "w").write("\n".join(out))
print("wrote", out_path)

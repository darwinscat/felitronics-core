# SPDX-License-Identifier: AGPL-3.0-or-later
# Dev-only: nulls what dump.cpp wrote against an independent numpy / scipy recompute. See README.md.
#
#   python check.py <dir with the dumps> <path to Mp3Window.h>
import re
import sys

import numpy as np
from scipy import fft as sfft
from scipy import signal
from scipy.ndimage import median_filter

d, header = sys.argv[1], sys.argv[2]
D = np.array([float(v) for v in re.findall(r'-?\d+\.\d{9}', open(header).read())])
assert len(D) == 512
GROUPS, FLOOR, CEIL, ZERO = 8, 1e-4, 1e4, 1e-2
failures = 0


def verdict(name, ok, text):
    global failures
    failures += 0 if ok else 1
    print('%s %s: %s' % ('PASS' if ok else 'FAIL', name, text))


# ---------------------------------------------------------------- the transforms, straight from their definitions
C = D / 32.0
POLY = np.cos((2.0 * np.arange(32)[:, None] + 1.0) * (np.arange(64)[None, :] - 16.0) * np.pi / 64.0)
MDCT36 = np.sin(np.pi / 36.0 * (np.arange(36) + 0.5))[:, None] * np.cos(
    np.pi / 72.0 * (2.0 * np.arange(36)[:, None] + 1.0 + 18.0) * (2.0 * np.arange(18)[None, :] + 1.0))
CI = np.array([-0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037])
CS, CA = 1.0 / np.sqrt(1.0 + CI * CI), CI / np.sqrt(1.0 + CI * CI)


def mp3_frames(x, offset):
    """x [n] -> [granules, 576]: the hybrid lines of every granule at sample offset `offset`."""
    phase, start = offset % 32, offset // 32
    count = (len(x) - phase - 512) // 32 + 1
    frames = np.lib.stride_tricks.sliding_window_view(x[phase:], 512)[::32][:count]
    sub = (frames[:, ::-1] * C).reshape(count, 8, 64).sum(axis=1) @ POLY.T          # [count, 32]
    g = (count - start) // 18 - 1
    blocks = sub[start: start + 18 * (g + 1)].reshape(g + 1, 18, 32).copy()
    blocks[:, 1::2, 1::2] *= -1.0
    win = np.concatenate([blocks[:-1], blocks[1:]], axis=1)                         # [g, 36, 32]
    xr = np.einsum('gts,tk->gsk', win, MDCT36).reshape(g, 576)
    out = xr.copy()
    for sb in range(1, 32):
        lo, hi = 18 * sb - 1 - np.arange(8), 18 * sb + np.arange(8)
        out[:, lo] = xr[:, lo] * CS + xr[:, hi] * CA
        out[:, hi] = xr[:, hi] * CS - xr[:, lo] * CA
    return out


def kbd(n, alpha):
    k = np.kaiser(n // 2 + 1, np.pi * alpha)
    c = np.cumsum(k)
    w = np.sqrt(c[:n // 2] / c[n // 2])
    return np.concatenate([w, w[::-1]])


def celt_window(frame=1920, overlap=120):
    w = np.zeros(frame)
    edge = (frame // 2 - overlap) // 2
    i = np.arange(overlap)
    rise = np.sin(0.5 * np.pi * np.sin(0.5 * np.pi * (i + 0.5) / overlap) ** 2)
    w[edge: edge + overlap] = rise
    w[edge + overlap: frame - edge - overlap] = 1.0
    w[frame - edge - overlap: frame - edge] = rise[::-1]
    return w


def mdct_frames(x, offset, window, count):
    """The MDCT by its definition through a DCT-IV of the folded frame (scipy's is twice ours)."""
    n = len(window)
    hop, q = n // 2, n // 4
    f = np.lib.stride_tricks.sliding_window_view(x[offset:], n)[::hop][:count] * window
    a, b, c, dd = f[:, :q], f[:, q:hop], f[:, hop:hop + q], f[:, hop + q:]
    u = np.concatenate([-c[:, ::-1] - dd, a - b[:, ::-1]], axis=1)
    return 0.5 * sfft.dct(u, type=4, axis=1)


KIND = {
    'mp3': (576, 576, None),
    'aac_sine': (1024, 1024, np.sin(np.pi / 2048 * (np.arange(2048) + 0.5))),
    'aac_kbd': (1024, 1024, kbd(2048, 4.0)),
    'celt': (960, 800, celt_window()),
}


def frames_at(pcm, kind, offset):
    hop, used, window = KIND[kind]
    if kind == 'mp3':
        return [mp3_frames(pcm[:, c], offset) for c in (0, 1)]
    x = pcm
    if kind == 'celt':
        x = pcm.copy()
        x[1:] -= np.float32(0.85000610) * pcm[:-1]
    count = (len(x) - 2 * hop) // hop
    return [mdct_frames(x[:, c].astype(np.float64), offset, window, count)[:, :used] for c in (0, 1)]


def signals(frames):
    l, r = frames
    return np.stack([l, r, 0.5 * (l + r), 0.5 * (l - r)])


def curve_of(pcm, kind):
    hop, used, _ = KIND[kind]
    out = np.zeros((hop, 4, GROUPS))
    ref = None
    for o in range(hop):
        sig = signals(frames_at(pcm, kind, o))
        if ref is None:
            ref = np.sqrt(np.mean(sig * sig, axis=1)) + 1e-30
        r = np.clip(np.abs(sig) / ref[:, None, :], FLOOR, CEIL)
        out[o] = (20.0 * np.log10(r)).mean(axis=1).reshape(4, GROUPS, -1).mean(axis=2)
    return out, ref


def summarise(curves):
    base = median_filter(curves, size=(17, 1, 1), mode='mirror')
    dip = (base - curves).astype(np.float64)
    local = dip.sum(axis=(1, 2))
    sigma = 1.4826 * np.median(np.abs(local - np.median(local))) + 1e-9
    o = int(np.argmax(local))
    n = len(local)
    far = np.ones(n, dtype=bool)
    for k in range(-2, 3):
        far[(o + k) % n] = False
    cell_sigma = 1.4826 * np.median(np.abs(dip - np.median(dip, axis=0, keepdims=True)), axis=0) + 1e-6
    return o, local[o] / sigma, local[far].max() / sigma, local[o], int(np.sum(dip[o] / cell_sigma > 4.0))


said = {l.split()[0]: l.split()[1:] for l in open(d + '/readings.txt')}
for kind in ('mp3', 'aac_sine', 'aac_kbd', 'celt'):
    hop, used, _ = KIND[kind]
    pcm = np.fromfile('%s/%s.pcm.f32' % (d, kind), dtype='<f4').reshape(-1, 2)
    theirs = np.fromfile('%s/%s.curve.f32' % (d, kind), dtype='<f4').reshape(hop, 4, GROUPS).astype(np.float64)
    ours, ref = curve_of(pcm.astype(np.float64) if kind != 'celt' else pcm, kind)
    diff = np.abs(ours - theirs).max()
    verdict('%s curve' % kind, diff < 2e-3, 'max |difference| %.2e dB over %d x 32 cells' % (diff, hop))
    coded_at, o, score, second, dipsum, cells = int(said[kind][0]), int(said[kind][1]), *map(float, said[kind][2:5]), int(said[kind][5])
    # the reading of THEIR curve, by numpy: the statistic alone, no transform in the way
    no, nscore, nsecond, ndip, ncells = summarise(theirs.astype(np.float32).reshape(hop, 4, GROUPS))
    ok = (no == o == coded_at and abs(nscore - score) < 1e-3 * score and abs(nsecond - second) < 1e-3 * max(second, 1.0)
          and abs(ndip - dipsum) < 1e-3 and ncells == cells)
    verdict('%s reading' % kind, ok, 'offset %d (coded at %d), score %.3f vs %.3f, second %.3f vs %.3f, dip %.3f vs %.3f, cells %d vs %d' % (
        o, coded_at, score, nscore, second, nsecond, dipsum, ndip, cells, ncells))
    sig = signals(frames_at(pcm.astype(np.float64) if kind != 'celt' else pcm, kind, o))
    zero = np.abs(sig) < ref[:, None, :] * ZERO
    share = zero.mean(axis=1).reshape(4, 32, -1).mean(axis=2).max(axis=0)
    zs = np.fromfile('%s/%s.zeros.f32' % (d, kind), dtype='<f4')
    verdict('%s zero share' % kind, np.abs(share - zs).max() < 5e-3, 'max |difference| %.2e over 32 bands, mean share %.3f' % (np.abs(share - zs).max(), share.mean()))

x = np.fromfile(d + '/original_left.f32', dtype='<f4').astype(np.float64)
for up, down in ((160, 147), (147, 160)):
    y = np.fromfile('%s/resample_%d_%d.f32' % (d, up, down), dtype='<f4')
    want = signal.resample_poly(x, up, down, window=('kaiser', 14.0))
    n = min(len(y), len(want))
    diff = np.abs(y[:n] - want[:n]).max()
    verdict('resample %d:%d' % (up, down), len(y) == len(want) and diff < 5e-6, 'length %d vs %d, max |difference| %.2e (peak %.3f)' % (len(y), len(want), diff, np.abs(want).max()))

print('ALL NULLS PASS' if failures == 0 else '%d NULLS FAILED' % failures)
sys.exit(1 if failures else 0)

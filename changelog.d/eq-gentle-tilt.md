### eq — a gentle tilt on first-order shelves

`matched::lowShelf1` / `highShelf1` (and their `…Db` forms) are first-order shelves, 6 dB/oct, beside `lowpass1` /
`highpass1`: the analog shelf with exactly half the gain in dB at f0, matched in magnitude at DC, at f0 and at
Nyquist, the pole in closed form. A `FilterType::Tilt` band with `slope == 6` now builds its tilt from these two
instead of the 2-pole shelves: the same pivot (unity at f0), the same ends (lows −gain, highs +gain), but the change
spreads over decades instead of the octaves around f0. At +3 dB about 1 kHz (48 kHz): 125 Hz −2.91 · 250 −2.64 ·
500 −1.79 · 1k 0 · 2k +1.79 · 4k +2.65 · 8k +2.91, where the 2-pole tilt reads −3.00 · −2.98 · −2.64 · 0 · +2.64 ·
+2.98 · +3.00.

Every other slope, the default 12 included, keeps the 2-pole tilt bit for bit. A Tilt band that already carries
slope 6 — say from a type switch in a product that keeps the slope across types — turns gentle with this release.

### io — a sample rate that rounds to 0 Hz is refused on both sides

The writer accepted a rate in (0, 0.5), rounded it to a 0 Hz header and wrote the file; the reader then returned
`sr = 0` with `ok = true`, and every division by the rate downstream became an infinity. The writer now refuses a
rate that rounds to 0, and the reader refuses a 0 Hz header with an error that names it. 0.5 Hz rounds to 1 Hz and
is written as before; no real audio rate changes.

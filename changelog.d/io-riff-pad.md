### io — the WAV writer pads an odd `data` chunk, so every file it writes is word-aligned RIFF

RIFF pads every chunk to an even length with one zero byte that the chunk's own size does not count and the RIFF
size does. `writeWav`/`writeWavMemory` never wrote that byte, so a 24-bit data chunk of an odd length — an odd
channel count times an odd frame count, 3 bytes a frame in mono — ended the file on an odd byte. Measured before the
fix: 24-bit mono, 1 frame → a 47-byte file (RIFF 39, data 3); 1001 frames → 3047 bytes (RIFF 3039, data 3003).
Forgiving readers took those as they were — ffprobe, afinfo, Python's `wave` and our own reader all read the audio —
but a chunk appended after one (a LIST/INFO tag, `bext`, iXML) starts on an odd offset. On such a file ffprobe read
the audio and silently lost the tag, and `readWav` rejected it as corrupt: a reader that skips the pad lands one
byte into the next chunk's header.

The writer now emits the pad after an odd data chunk and counts it in the RIFF size; the data size stays
frames · blockAlign. What a caller sees change: only a 24-bit image with an odd channel count AND an odd frame count
differs — one zero byte longer, with a RIFF size one larger — so code that took the file size to be 44 + data must
allow for it. Every 16-bit and 32-bit float image, and every other 24-bit one, is byte-identical to before. The
reader was already right — it steps over a pad after any odd chunk, and still accepts an odd final chunk whose pad
is missing, which is how every file written before this one reads — and is unchanged but for a comment.

New suite `felitronics_io_riff_tests` (35 checks). 16/24-bit PCM and 32-bit float × 1/2/3 channels × 1/2/3/1000/1001
frames: an even file, RIFF = file − 8, data = frames · blockAlign, a zero pad exactly when the data is odd, a strict
chunk walk that covers the file exactly, write∘read and read∘write identities, and a LIST appended after the image
reading back with the audio intact. The reader over hand-built images: odd LIST/ID3/junk chunks before `fmt `,
between and after an odd data chunk, and the old unpadded layout. Planted failures run through the same
instruments: the old layout fails the walk and the append on exactly the six odd cases and passes the other 39, and
a walk that does not skip the pad loses the four images with an odd chunk before another. Against the old `Wav.h`
the suite reads 7 of 35 red. `felitronics_io_grid_tests` now takes the data length from the header and requires the
pad, so its 24-bit image (4869 codes, 14607 bytes of data) exercises it too; its planted 24-bit index checks no
longer read past an image that failed to decode, which crashed the suite instead of failing a check. Cross-checked
outside the suite: the 45 images and their LIST-tagged twins read correctly in ffprobe, ffmpeg, afinfo and Python's
`wave`, the tag seen on every one.

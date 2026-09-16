<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools · wasm — loudness and clipped runs of a stream, read between the pieces

**`fc_probe` gains a streaming surface** for two instruments that were already streaming:
`analysis::DeterministicLoudnessMeter` and `analysis::ClipDetector`, behind one class, `fcore::StreamProbe`
(`tools/fcore_stream.h`). A page opens a stream with `fc_stream_create (rate, channels)`, feeds planar PCM
with `fc_stream_process`, and between pieces reads `fc_stream_loudness` (momentary, short-term and integrated
LUFS, samples consumed, gating blocks dropped past the one-hour store), `fc_stream_clips_count` and
`fc_stream_clips (h, from, out, cap)` — the runs decided so far, polled by index. `fc_stream_finish` decides the
last runs; `fc_stream_destroy` frees the stream. Nothing in either instrument changed.

**Handles, not a singleton**, so two streams run at once — at most 16. A handle is a serial looked up in a table
and never reused, so a stale, forged or failed handle is refused rather than dereferenced. Buffers are checked
as the rest of `fc_probe` checks them, `cap` is in doubles and the return in runs as in
`fc_probe_clips_runs`, and a refused piece that carried samples poisons the stream: every reader answers 0.

**The deterministic meter, not the system one** that `fc_probe_run` uses: the system meter's block energies
differ between native and wasm at 8000, 88200 and 192000 Hz. `fcore_measure stream` prints the same bytes as
`tools/wasm/stream-parity.mjs`, and CI diffs them at 44100, 48000 and 88200 Hz on the release and checked
modules. `felitronics_stream_abi_tests` feeds pieces of 1, 4096 and random sizes through three handles at once
and requires the same integrated loudness, bit for bit, and the same runs as the same instruments given the
whole buffer in one call.

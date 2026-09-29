### limiter — the K13 peak clipper's reduction, sample by sample

`TruePeakLimiterTap` gains `peakClipReductionDb`: one positive dB value per OVERSAMPLED sample, on the same grid as
`gainReductionDb` and `linkedPeakLin`, zero where the clipper did nothing. Until now the clipper published only
aggregates — the deepest reduction, the run count, the histogram — which say how much it clipped but not where; a
mastering report that draws the clip on the programme's time axis needs the trace. The field is appended, so every
existing `{ gainReductionDb, linkedPeakLin, capacity }` initialiser still compiles; a short `capacity` refuses the
call with this pointer alone, as it does with the other two. Tapped and untapped calls keep the same PCM and the
same aggregates, bit for bit.

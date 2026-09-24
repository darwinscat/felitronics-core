### Fixed — a second run of `StereoBandBursts` lost the other axis's witness

`reset()` cleared every list and every latch except the **commit cursor**. After a run `committed_` equals
the event count, so on the next run `while (committed_ < stored)` fires for none of the first N events and
their cross readings stay at their default: power 0, eligible false, hop −1 — which this class's own header
calls *"reads as 'the other axis was silent there'"*. The first answer was right and the second was
confidently wrong.

A page holding one probe per request never saw it. A witness run on a render's **output**, and every test
that reuses an instance, did — which is how the consumer found it on v0.48.0, reproduced twice before
reporting. The core's own suite missed it for one reason: it never ran a single instance twice.

The regression fixture runs the same object three times and demands identical bits, on a **hard-panned**
programme on purpose — with content only in L, Side equals Mid exactly, so the cross reading equals the
event's own `peakPower` and an unwritten one is off by the whole value rather than by a little. Three runs
and not two, so the cursor has to be cleared rather than merely decremented. Planting the shipped behaviour
reproduces it exactly: cross **0.000000** on the second run against **0.010104** on the first.

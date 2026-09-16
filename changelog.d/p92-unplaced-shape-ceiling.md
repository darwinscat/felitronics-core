<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam — what the ledger cannot place costs ONE allowance: never zero, never the face value

`felitronics::nam::detail`'s receptive-field registry promises an **upper bound** on a capture's memory,
and two readers spend its answer: law 11a's drain and P47's stream restart feed a lane that much digital
silence. Where it could not read a shape it answered **zero**, which is the tail of the previous sound
coming out of silence. It now answers what it read, plus one named allowance.

- **The measured case is NAM's own shipped capture.** `example_models/slimmable_wavenet.nam` with its
  config moved under `config.model`, and a top-level `layers` carrying nothing but the slimmable marker,
  **loads** (NAM's WaveNet parser delegates on the marker and the delegate reads `config.model`), makes
  the same sound, and its impulse reaches sample 2046 either way. The registry answered **2047** for the
  flat file and **0** for the wrapped one, and NAM's own answer for the architecture is `return 0`. All
  three questions were blind together: an LSTM conditioner inside the wrapper read as not recurrent, and a
  dense `Linear` one lost its 2048-sample ring.
- **Three events mean "cannot place", all read off the config and none off NAM's dispatch:** an object
  carrying a model's vocabulary under a key the registry does not read (at a config, a layer entry or a
  submodel entry, singly or in an array — under any key name); a value that is there and cannot be read
  (a slimmable dilation of `4294967396`, which NAM builds as 100, answered 0); and a reading the registry
  sets aside (a declared field beside a stack). Each adds `kUnreadShapeCeiling` = **48 000 samples**,
  **once per tree**, to what the registry did read, and charges the ring. It fires on **none** of the
  1229 distinct captures on the author's machine; the flat shipped slimmable still reads 2047.
- **Why once, why added, why not the face value** — each is a measured failure of the alternative. An
  allowance per node turned a 1.6 MB file of 30-byte dead siblings into an **INT_MAX** drain; a max at the
  root let a known 100 001-sample stack swallow the allowance of an unreadable stage in series with it;
  and trusting the face value of a dead `"receptive_field": 2147483647` beside a real Standard's stack —
  a file NAM loads unchanged — made `reset()` run for **about half an hour per lane** (2³¹ samples at 0.79–0.89 µs each).
- **The price, measured through `NamStage::reset()`** on real captures rewrapped, per lane at 48 kHz:
  **39.7 ms** on `wavenet_a1_standard` at a 256 block (44.5 ms at 64), **24.1 ms** on A2's submodel,
  **4.5 ms** on the slimmable itself. On shapes that ship: identity.
- **What moved, on purpose, and only in synthetic rows:** a stale `receptive_field` beside a stack
  (9 → 48 009), a stray top-level `dilations` beside `layers` (2 → 48 002), refused values (1e300,
  `4294967396`, a string) from "absent" to the allowance. **Closed on the way:** a 5000-tap `Linear`
  carrying a readable stray `layers` array answered 2 (reach 4999; now 48 002), and a `ConvNet` carrying a
  non-empty dead `layers` array answered 0 (NAM: 256; now 256).
- **What it does not close, stated rather than hidden — two doors, and shutting either opens the other.**
  A LIVE memory the registry cannot place and that is LONGER than the allowance drains short by the
  difference (a wrapped model with an inner field of 60 000; a 60 001-tap `Linear` whose declared field is
  set aside beside a dead `layers` array) — the longest real field is 6347. And a DEAD number the registry
  PLACES is still trusted at face value, exactly as before: a lower reading with no stack beside it (a dead
  `dilations:[2e9]` beside a `Linear`'s declared field), the wrapped form's own decoy stack, a `layers`
  array on an architecture that never reads it. Which door stays open is registered as a policy question.
  A differential fuzz of 2 000 loadable configs against NAM found no short answer outside these two and
  the recurrent exception. Separately, NAM's own recursive copy of the config takes the host down on a
  file about 2 000 levels deep before the registry runs. The only new recursion — looking for an
  architecture through unplaced nodes — is bounded at 32 unplaced hops; the nesting the registry already
  walked is read exactly as before, and a guard on it was measured to drain a 65-deep conditioner chain
  short.

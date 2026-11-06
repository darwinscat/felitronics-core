### basscomp — a new module: `BassComp`, the bass amp's COMP

`felitronics::basscomp::BassComp` moved here from OrbitAmp (`src/core/BassComp.h`), JUCE-free and its arithmetic
unchanged — header-only, depends on `core`, `eq` and `dynamics`; a module of its own because it is composition
(`dynamics` stays independent of `eq`). EACH CONTROL ONE JOB: the threshold (-48 … 0 dB) and the ratio (1.5 … 10:1)
are their own; the knob is BLEND, parallel compression — 0 is no compressor at all (once glided there nothing is
computed and the block passes bit for bit), 1 all of the signal compressed, sample-aligned (no lookahead). Make-up is
automatic and level-neutral: what the squeeze has taken from this playing, averaged over 400 ms while a note sounds
and held through pauses. Five characters: core's compressor on the RMS (SMOOTH) or the peaks (PUNCHY), and
`dynamics::CharacterComp`'s OPTO, FET and SQUASH, a switch between the two compressors crossing over in 15 ms. The
detector hears the bass through a 30 Hz, 12 dB/oct high pass (`eq::CutBiquad`). Every move glides (20 ms for each
control's whole travel), and the controls reach the compressor on a 1 ms cadence counted from the prepare: the same
playing gives the same samples whatever the host's blocks. The standby is `core::BypassFade`'s span, handed in per
block: standing by with BLEND up it is fed, not heard.

What the plug-in used to write into its picture from inside `process` is now a read-out: with `pictured` set, a
heard block leaves a `Picture` (the signal handed in, the squeezed signal before the make-up, the same blended as
the sound is, the squeeze per sample, the detector, whether the trace breaks before it) until the next call;
`picture()` is null where there is none. Nothing of it reaches the sound.

Measured (`modules/basscomp/tests/BassCompTests.cpp`): BLEND 0 is the note bit for bit; a steady note at the
reference -18 dBFS RMS keeps its level within 1.2 dB in all five characters at BLEND ½ and 1, and at -40 dB, 8:1; a
-40 dB note under a -12 dB threshold stays within 1e-5 of itself; host blocks of 512, 64 and 333 play the same
samples (SMOOTH, OPTO, FET); standing by and settled, the dry bit for bit; a prepare's `snapAmount` is heard from its
first block and a recall up from 0 is fed, not heard, through its warm; the picture is there only while pictured and
heard and says when its trace breaks; the indicator reads the squeeze heard; no allocation in 40 blocks with a
character switch among them.

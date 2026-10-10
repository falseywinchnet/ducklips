# Smooth jet-separation diagnostic 037 — 2026-10-10

Change only the hard attached/separated force switch to a smoothstep across
local area/minimum area ratios 1.0–1.4, centered on the previous 1.2 threshold.
The papa/mama design and speech timing are identical to baseline 036.

Fresh raw audio: `out/words-papa-mama-37/outlet.wav`.
SHA-256: cb95f5489ca65525d7a75e66f5e1f6ee905422e28a87590aa91f6fd3a9594868.
Source snapshot and full enhanced analysis accompany that take.

The native-rate cycle diagnostic in `out/word-cycle-baseline-smooth.png`
compares all eight axial cord heights, valve flow, approximate applied wall
pressure (body pressure plus its offset), and raw output during 40 ms of the
first vowel. Both versions show irregular closures and large pressure
excursions. Smoothing the separation transition does not remove those
excursions. Do not promote this as a voice-quality fix. No listener preference
for this take has been established.

The force diagnostic is a snapshot of two staggered quantities, not an exact
energy audit. In particular the original below-mouth pressure calculation
and the added nonlinear valve both influence the fold; the offset subtracts
a low-pass mean, leaving the original acoustic pressure present. The evidence
points toward investigating contact and pressure coupling before more tuning,
but does not isolate a unique cause or prove a repair.

The recognized default remains the hard-switch baseline. It is protected by
an exact raw-WAV hash test. `node web/tube/render-speech.mjs 'papa mama'`
uses that baseline. A final `smooth-separation` argument enables this experiment
when a design path is also supplied.

The complete 428 by 684 native-cell magma enhanced field was inspected. It
retains four broad vowel groups and irregular within-group energy; it does
not demonstrate a clean harmonic source. No grayscale or ordinary STFT
comparison was used.

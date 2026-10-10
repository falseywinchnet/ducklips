# Flat-tube restart — 2026-10-10

The user rejected the preceding five Papa experiments and requested a fresh
geometry-driven tube engine and a directly editable 2D web simulation. They
clarified that the editing view is a side cutaway with two walls, not an
unrolled circumferential surface. This implementation follows that view.

## Scope delivered

Independent browser synthesis core; six bands, four nodes per surface per
band; support-derived material modes; upper/lower support dragging; additive
volume blocks with normalized enclosed volume; mass/elasticity/material
editing; bidirectional axial/transverse/crisscross air links; separate solid
links with material-dependent FIR transmission; slight reservoir-pressure
wander; local pressure/wall coupling; outlet sampling and real Web Audio
monitoring. It includes undo, save/open, live WAV recording and automatic
complete enhanced depiction through a localhost-only server.

The new core imports no previous Ducklips synthesis code. It contains no
prescribed voice waveform, phoneme timeline or recording playback. The existing
enhanced renderer is reused for observation. Earlier work is preserved at
archive/pre-flat-tube-20261010 (36d670d), with README.legacy.md identifying the
historical entry point. The prior five hypotheses are now recorded as rejected.

## Visual rounds

All paths below are local generated artifacts under out/, intentionally not
committed. Every listed audio take has its own complete enhanced field.

1. `out/flat-tube-first`: new tube, 0.8 s. Wide rough activity and large late
   bursts. Peak 4.983 before monitor protection. Rejected as a transport defect.
2. `out/flat-tube-normalized`: delay states changed to energy-normalized pressure
   waves so an edited/changing admittance does not invent stored wave energy.
   The full image has organized harmonic ridges instead of the former large
   irregular columns. This is a numerical representation correction, not a
   successful voice judgment.
3. `out/flat-tube-current`: same corrected physical core, fused propagation FIR,
   compact delay storage and fixed outlet conversion gain. Fresh 0.8 s, 48 kHz
   WAV from a 24 kHz engine. Peak 0.12667; RMS 0.01399; no model fault. Enhanced
   field inspected in full (684 x 151 computed cells, hop 256, original magma).
4. `out/flat-tube-web/2026-10-10T18-05-06.594Z-5adfc0`: first actual browser
   capture, 2 s at 48 kHz, raw peak 0.16925 and RMS 0.05509. Full 684 x 376 field
   inspected at one pixel per cell. Strong stable harmonic ridges, with slowly
   varying positions and levels; a sustained periodic instrument output.
   Human naturalness is unassessed. This UI-path test predates automatic
   design/source receipts; do not infer a frozen source receipt for this take.
5. `out/flat-tube-web/2026-10-10T18-07-43.484Z-e7e49a`: final browser recording
   path with raw WAV hash, frozen core/worklet/UI source, initial design and
   any edits during capture in take.json. See validation.json for the final
   visual review and checks.

There was no listening-tool audition. Claims above are visual or measured,
not claims that the model sounds human. The monitor was set to zero during
browser automation; raw recordings bypass that monitor control.

## Verification and limits

Ten focused tests cover block volume conservation, eigenpair residuals,
mass/elasticity dependence, anchor/volume-induced mode changes, silence at
unpressurized rest, causal/bidirectional connectivity, physical-state projection
when supports change, deterministic replay, invalid parameters and ringdown
after reservoir removal. A ten-second default run remained finite; a subsequent
three-second release reduced late output RMS to 3.43e-8. These establish the
checked implementation properties, not physiological correctness.

Browser interactions verified volume insertion, upper-anchor dragging (opening
0.54 -> 0.35 mm with modal changes), keyboard adjustment, material switching,
undo, JSON save/open, live audio processing, fresh capture and the enhanced
request. Browser error/warning log was empty during the final capture.

Initial core speed fell below 10x. Fractional delay interpolation and attenuation
were combined into one four-tap FIR, and fixed overlarge history buffers were
replaced by path-sized buffers. No material modes were removed. Final Node
core runs were 11.38x, 10.15x and 10.38x; the browser showed 10.1–10.8x. These
measure the 24 kHz physical core, excluding geometry compilation, GUI rendering,
48 kHz reconstruction/monitoring and offline enhanced analysis. See benchmark.json.
The 0.8 s full enhanced render took about 23 s; the first 2 s live field took
123.5 s. The spectral tool's cost is explicitly outside the realtime engine.

The model still needs physical validation of moving geometry, propagation
speed/mesh dispersion, contact, and glottal flow. It does not resolve the
separated jet or a validated viscous constriction law. Complete airtight
closure is approximated by a small finite opening. Material parameters and
solid-wave impedance are prototypes. No word or natural voice is accepted.
These limits are documented in web/tube/MODEL.md rather than hidden behind
successful unit tests or a regular spectrogram.

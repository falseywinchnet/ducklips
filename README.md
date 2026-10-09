# ducklips

A first-principles **human speech** engine. Stream articulatory and physiological
nudges into a body; its moving airway, self-oscillating vocal folds, airflow and
acoustic radiation produce the sound. Body geometry and tissue properties are
configuration, not a speaker embedding or a bank of recorded sounds.

**Status: active research, not perfected speech.** Current work connects a
volume-conserving tongue body to the canonical moving acoustic network and
fixes physical and pronunciation-planning failures. Local material transport
now resolves the tested Foxtrot attachment miss while retaining the dorsal
contact and lower anchor. Human listening of that alphabet found 14 words
legible but mechanical, six distorted, two with clicks and four unintelligible;
see the [listening record](reports/LISTENING_2026_10_09.md). The canonical
renderer now removes [unphysical history in fully closed passages](reports/CLOSED_PASSAGE_HISTORY.md).
The preceding listening result does not score subsequent renderings. The
[spatial pressure experiment](reports/SPATIAL_PRESSURE.md) now connects local
pressure storage and airflow to this engine. It reduces one Foxtrot transient
but worsens another, exposing a discontinuity in the midpoint airway-area
measurement when the tip passes a section. It remains opt-in research.
Performance optimization
is parked by user instruction. Numerical stability and a plausible enhanced
spectrogram do not establish intelligibility or human realism. See the
[current attachment model](reports/LOCAL_TONGUE_ATTACHMENT.md), its
[speech integration](reports/CONSERVATIVE_TONGUE_SPEECH.md), and the
explicit [research gates](RESEARCH.md). Earlier engines and their performance
measurements remain research references; they are not acoustically equivalent.

Listening follow-up on the spatial-pressure recording: the user reports
[November improved, with Foxtrot still jammed around /tr/](reports/LISTENING_SPATIAL_2026_10_09.md).
A new [solid-airspace geometry reference](reports/PLANAR_AIRSPACE.md) measures
air above and below the tongue, retaining the /t/ closure while eliminating
the tested finite sampling jumps. It is independently checked against Boolean
geometry and is not yet a speech-rendering path. The native
[mean-pressure solver now retains cyclic and disconnected air routes](reports/CAVITY_GRAPH.md);
fixed Foxtrot snapshots preserve the sealed /t/ chambers. Moving-state and
acoustic integration remain unfinished. [Zero-storage transport checks](reports/CAVITY_TRANSPORT.md)
now conserve cell contents through disappearance. [Separate regions within
cells](reports/PLANAR_REGIONS.md) resolve the split-cell contacts at all 66 saved
Foxtrot poses and preserve the /t/ seal in a pressure control. [Sampling moving
contours](reports/PLANAR_MOTION.md) now passes an independent 521-pose audit and
exposes openings invented by endpoint-metric interpolation. [Opening histories
between poses](reports/PLANAR_FACE_MOTION.md) now retain within-step contact and
resolve a false post-closure pressure transfer in a controlled fixture.
[Stored-state continuation](reports/REGION_STATE_CONTINUATION.md) now preserves
air behind a moving seal through graph reindexing and audits correspondence
across all 520 saved subintervals. [Partitions through the Foxtrot contact
points](reports/CONTACT_PRESSURE_PARTITION.md) now retain finite neck coupling
without an instantaneous pressure merge, and the native solver now checks
small-region conservation against each region's own scale. A native
[moving-pressure replay](reports/PLANAR_PRESSURE_MOTION.md) now combines actual
contour motion, opening histories and conserved storage through two separately
initialized closure/release intervals. It preserves the sealed chambers, but
the coarse release control reaches a 12% density reduction; its prescribed
mobility and pressure magnitudes are not physically validated. Complete
contact-event evolution, physical interval flow and acoustic integration
remain required for speech.

The subsequent [inertive airspace replay](reports/PLANAR_INERTIVE_MOTION.md)
now connects these moving volumes and openings to the canonical pressure-and-
momentum integrator. The opening's air inertia replaces the control's arbitrary
resistance. In the reopening fixture, the worst density reduction falls from
12% to 3.7%; the remaining shared front-mouth pressure drop is consistent with
its sealed lips and expanding volume. This is a coarse inviscid diagnostic,
not a new speech recording or a validated flow law for the human mouth.

[Explicit domain-cap attachments](reports/PLANAR_BOUNDARY_PORTS.md) now retain
each opening's owning region, area, centroid and outward normal, including
contact with exterior tissue. Connecting the mouth-end cap to an ideal ambient
pressure boundary reduces the reopening control's worst density reduction
from 3.7% to 0.65%, while retaining the internal tongue seal. Actual lip shaping,
vocal-fold coupling and radiation losses are still required; a domain cap is
not the complete anatomical opening.

The [general-graph acoustic reference](reports/MOVING_GRAPH.md) now retains
branch momentum through cycles, allowing multiple mouth regions to share one
outlet volume. Independent discrete and continuous controls check this
coupling; no new speech has been rendered with it.

The [articulated-lip reference](reports/PLANAR_LIP_COUPLING.md) now implements
those replacement seams and connects the final lip section areas and lengths
from the native articulator replay. Both front-mouth regions share the same
lip tube; its narrowing changes storage and drives airflow. The tested open-lip
intervals conserve stored air and retain the tongue seal. Closed-lip transport,
laryngeal coupling and radiation are still required for this to render speech.

The [lip contact investigation](reports/LIP_CONTACT_STORAGE.md) now separates
actual air volume from the contact aperture. In all 562 NATO frames with a
zero final-section midpoint area, the existing axial lip law still leaves
air in that section's inner shoulder. Native integration retains that air
while sealing the outlet, and coupled controls preserve its pressure through
contact.

The [coupled profile adapter](reports/PLANAR_LIP_PROFILE.md) now uses those
separate measurements in the native planar reference. Moving Bravo closure and reopening
controls carry mouth and lip state through successive intervals, retain the
air behind the seal, and conserve stored air including exterior flow. Those
initial controls use a sealed inlet; no new listening result is attributed
to them.

The [throat and source connection](reports/PLANAR_SOURCE.md) now adds the
existing epilarynx sections and an implicit held-aperture glottal inlet to
that same stored-air network. It includes source-path inertia and pressure
feedback, with independent mass, energy and continuous-equation checks.
The moving inviscid controls still need temporal convergence; this is not
yet a speech path.

[Moving vocal folds](reports/GRAPH_PHONATION.md) now couple to that source
with pressure feedback within each step and complete state continuation.
Independent calculations agree, but a work audit exposes an unresolved
exchange between the inherited fold forces and the algebraic glottis during
closure. Internal glottal storage/displacement and matching pressure work
need resolution before this source can be treated as physically complete.
No new listening result is attributed to it.

[Finite glottal storage](reports/GLOTTAL_STORAGE.md) now gives the moving
folds two air chambers whose compression, mass transfer and tissue forces
obey one work balance. Independent thermodynamic and continuous-motion
checks pass. Its flows are still supplied inputs; connecting passage
momentum, jet losses and the tract remains necessary before speech.

[Pressure-driven chamber coupling](reports/GRAPH_GLOTTAL_STORAGE.md) now
solves those three flows against retained passage momentum and the tract
load. Independent continuous controls pass, including reverse flow and
undriven decay. The passage coefficients remain fixed controls; their
dependence on moving fold geometry and full contact is still required.

[Mass- and geometry-derived passages](reports/GLOTTAL_PASSAGE_GEOMETRY.md)
now replace those fixed coefficients in an explicit articulated mode.
Kinetic head, changing inertia and viscous losses participate in the same
work balance. Constant reservoir pressure produces repeated fold motion in
a control; narrow-volume compression, contact and the loss model remain
physical questions, and long-run temporal accuracy needs finer stepping.

[Near-contact diagnosis and contour geometry](reports/GLOTTAL_CONTACT_GEOMETRY.md)
now show that strong compression persists under refinement while air escapes
from both ends. The new geometry keeps wedge-shaped storage behind an exact
seal and derives force areas and passage integrals from the same surface.
It is independently checked but not yet coupled into the source dynamics.

This repository preserves the supplied research snapshot as its first commit.
[DESIGN.md](DESIGN.md), [FINDINGS.md](FINDINGS.md), and
[surface/NOTES.md](surface/NOTES.md) preserve earlier mechanisms, measurements,
listening feedback, and rejected ideas. Their earlier machine/path references
and performance values are historical, not fresh validation of this checkout.

## Current speech integration

`ducklips_say --canonical-motion` and `ducklips_stream --canonical-motion`
select the new full canonical moving-network path. It connects branch-history
transport, work-conjugate glottal feedback, viscothermal memory and a passive
radiation load to the articulators and vocal folds. The old compact modal
path remains the default; both CLIs print their selected path. This reference
is slow, and pronunciation/realism remain unverified. Performance work is
deferred by user instruction while the integrated speech path is developed.
See [CANONICAL_SPEECH.md](reports/CANONICAL_SPEECH.md) for the implementation,
rejected instability and remaining physical limitations.

Add `--conservative-tongue` to select the new material body deformation and its
attached blade. Both flags are required. It reports failed geometry solves,
sampled body/wall overlap, and unmet tip travel separately. A failed geometry
solve stops rendering; it does not substitute the rigid body. The blade still
has a fixed length; local body deformation now advances its attachment when
needed. This is conservative kinematics, not calibrated muscle dynamics.
See [LOCAL_TONGUE_ATTACHMENT.md](reports/LOCAL_TONGUE_ATTACHMENT.md).

Add `--spatial-pressure` to either canonical CLI to exercise pressure storage
on the airway tree, local constriction flow, and separate mouth/nose outflow.
It requires `--canonical-motion`. Solver failure stops the render explicitly.
This experiment is not a new listening-approved voice: its remaining geometry
transients are documented in [SPATIAL_PRESSURE.md](reports/SPATIAL_PRESSURE.md).
For diagnostics, `ducklips_say --until SECONDS` renders a prefix of the full
plan without changing its preceding cues. The geometry replay tool accepts
`--airway CSV` to export every oral section's area and length at each block.

## Build

Requires a C++20 compiler and CMake 3.24+. The core has no downloaded dependencies.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The Python reference/analysis tools use NumPy, SciPy, Matplotlib, and SoundFile.
On this machine they are already available in
`/Users/joshuahkuttenkuler/miniforge3/bin/python`.

## Speak and inspect

```sh
./build/ducklips_say --canonical-motion --conservative-tongue --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Alfa. Bravo. Charlie.' --out out/nato.wav
./build/ducklips_say --canonical-motion --conservative-tongue --voice alto --rate 1.25 \
  --phonemes 'M AY1 K | S IY0 EH1 R AX0' --out out/phones.wav
python tools/enhanced_specgram/render.py out/nato.wav --out out/nato-enhanced --max-hz 8000
```

Create `out/` before writing there. Stock bodies: `baritone`, `alto`, `child`,
`gravel`, `whisper`. Explicit dictionary definitions override built-ins;
unknown words use a deliberately limited letter-to-sound fallback and are
reported. `--report` shows the phone plan and vowel fits. `--probe` and `--trace`
expose internal physics. `--normalize` is optional and is **off** for measurements.
The benchmark exits nonzero on any render below 10× realtime, nonfinite sample,
or clipped sample before output limiting. It reports prepare and total wall
time separately. The enhanced spectrogram is an offline diagnostic, outside
the audio rendering budget.

## Stream nudges directly

The text planner is optional. `ducklips_stream` reads commands incrementally
from stdin and writes **48 kHz stereo little-endian float32** to stdout. Logging
goes to stderr. `render UNTIL` emits through an absolute score time; subsequent
commands continue the same body and wave state.

```sh
./build/ducklips_stream --canonical-motion --conservative-tongue --body bodies/reference.body --speed 1 \
  < examples/diphthongs.nudges > out/diphthongs.f32
```

Example protocol:

```text
0 lungs 600 40
0 glottis 0.000003 170
0 vowel AA 75
0.2 diphthong AY 0.75 75
render 0.6
0.75 pitch 108 60
1.1 diphthong AW 0.75 75
2.0 lungs 0 40
render 2.5
```

Commands use `TIME KIND ... RATE`. Times/durations are seconds, pressure is Pa,
areas are m², pitch is Hz, and motion rates are rad/s. `--speed` scales score time
and requested motion; tissue speed limits still apply. `speaking_rate` in a body
file governs the text planner; the raw nudge protocol uses its explicit speed.
Diphthongs specify onset and offglide targets; the articulator dynamics generate
the motion. The current offglide starts at 55% of the requested duration, an
explicit timing convention, not a measured universal physiological law.

Other commands: `constrict PHONE ACTIVATION RATE`, `release lips|tip|body RATE`,
and `pitch`, `lungs`, `glottis`, `velum`, `effort`, `loudness` followed by value
and rate. Consonant nudges change constriction geometry; the host supplies
voicing and breathing separately. There is no automatic phrase script in this
interface. Stale cues, bad values, and queue overflow are rejected. The control
resolution is 1 ms, and cues within a block are applied at that block's start.

The C++ adapter is [NudgeStream](include/ducklips/nudges.hpp). `submit` and
`render` require one serialized owner; rendering allocates no heap storage.
For a separate producer and audio thread, use the underlying
`Synthesizer::enqueue` SPSC interface with chronologically ordered `Cue` values.

## Bodies and references

[Body files](bodies/reference.body) are strict `key=value` configurations in SI
units. They expose tract and nasal geometry, fold masses/stiffness/dimensions,
pressure, pitch habits, physiological variation, and individual articulator
speeds. Unknown fields, duplicates, invalid numbers, and inconsistent glottal
areas fail transactionally. Bounds guard this model's useful domain; they are
not asserted anatomical population limits. A body file is not evidence of a
successful individual-speaker match. [body_config.cpp](src/body_config.cpp)
contains the supported fields and ranges.

```sh
python tools/fetch_samples.py --prepare
```

The [reference corpus](samples/README.md) is human speech only: a hash-pinned
NATO alphabet/digits reading, plus 26 attributed word excerpts. Downloads are
bounded; raw and derived audio stay outside Git. There are no model weights,
recorded synthesis units, singing, music, or nonhuman vocal references.

## Enhanced speech field

The preserved full double-inverse / registered-aperture pipeline is in
[tools/enhanced_specgram](tools/enhanced_specgram/README.md):

```sh
python tools/enhanced_specgram/render.py out/nato.wav \
  --max-hz 8000 --channel 0 --out out/nato_field
```

It saves PNG, NPZ, and a JSON receipt with exact axes, parameters, source hash,
and timing. The PNG shows the complete enhanced field alone, using the original
magma colors and display range. Every computed row is preserved. The enhanced
field is not an invertible spectrum, calibrated energy measure, or realism score.

Inspect fine structure with the native-cell viewer, not the overview PNGs:

```sh
python tools/native_fields.py \
  out/enhanced_specgram/nato_human.npz \
  out/enhanced_specgram/nato_generated.npz \
  --out out/enhanced_specgram/native-fields.html \
  --png-prefix out/enhanced_specgram/full-enhanced
```

This emits an inline visualization fragment plus optional complete color PNGs.
The viewer displays every supplied row and time column of the enhanced field,
with at least one CSS pixel per row and column. The plot grows vertically when
rows are enlarged. There is no cropped viewport and no ordinary-STFT panel.
The original saved amplitude reference and dB range determine the original
magma colors; lossless PNG embedding preserves those colors exactly. The NPZ
retains full-precision linear fields and the original axes. Ordinary STFT and
intermediate fields remain available for targeted numerical investigations,
without being added to the inspection view by default.

`tools/compare_fields.py` is only a coarse timing overview. Setting
`interpolation='nearest'` does **not** preserve details when a dense field is
drawn in a smaller panel or the resulting PNG is reduced again for display.
Never use that overview to judge the enhanced field's additional fine detail.

## Detailed surface experiment

```sh
cmake -S surface -B build/surface -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/surface
./build/surface/surface_say --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Alfa.' --out out/surface.wav
```

This is the slow research reference, not a 10× rendering solution. `--tube`
removes transverse surface structure for diagnosis; it is not an accepted
fidelity-preserving replacement. The optional GUI.Forms adapter requires an
external checkout passed through `GUI_FORMS_SOURCE`; it was not built in this
standalone validation. `boundary/` preserves the earlier boundary experiments.

[Frozen port-response probes](surface/FROZEN_PORTS.md) provide an independently
verified frequency-domain equivalent of the detailed network for developing
and checking a faster physical representation. It includes transverse and
floor paths; it does not yet replace the speech renderer.

Third-party enhanced-field source keeps its original license and provenance.
Human reference audio and segmentation have separate licenses and attribution
under `samples/`. No new blanket license has been imposed on the supplied code.

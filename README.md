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
now conserve cell contents through disappearance, while a full 66-pose audit
exposes split cells and unresolved pressure refinement that block speech integration.

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

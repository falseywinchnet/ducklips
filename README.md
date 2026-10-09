# ducklips

A first-principles **human speech** engine. Stream articulatory and physiological
nudges into a body; its moving airway, self-oscillating vocal folds, airflow and
acoustic radiation produce the sound. Body geometry and tissue properties are
configuration, not a speaker embedding or a bank of recorded sounds.

**Status: active research, not perfected speech.** The compact modal engine
passes the local 10× realtime benchmark; the detailed surface engine does not.
The engines are not acoustically equivalent. The surface model has a mouth
floor and transverse propagation that the compact engine does not yet retain.
Neither numerical stability nor a plausible spectrogram proves intelligibility
or human realism. See [current findings](reports/2026-10-09.md) and the explicit
[research gates](RESEARCH.md).

This repository preserves the supplied research snapshot as its first commit.
[DESIGN.md](DESIGN.md), [FINDINGS.md](FINDINGS.md), and
[surface/NOTES.md](surface/NOTES.md) preserve earlier mechanisms, measurements,
listening feedback, and rejected ideas. Their earlier machine/path references
and performance values are historical, not fresh validation of this checkout.

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
./build/ducklips_say --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Alfa. Bravo. Charlie.' --out out/nato.wav
./build/ducklips_say --voice alto --rate 1.25 \
  --phonemes 'M AY1 K | S IY0 EH1 R AX0' --out out/phones.wav
./build/ducklips_analyze out/nato.wav --spectrogram out/nato.png --clicks
python tools/benchmark.py
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
./build/ducklips_stream --body bodies/reference.body --speed 1 \
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
  --start 0 --duration 1.5 --max-hz 8000 --channel 0 --out out/nato_field
```

It saves PNG, NPZ, and a JSON receipt with exact axes, parameters, source hash,
and timing. Keep its ordinary STFT beside the enhanced field: registration and
reassignment improve inspection but do not make an invertible spectrum, a
calibrated energy measure, or a realism score.

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

Third-party enhanced-field source keeps its original license and provenance.
Human reference audio and segmentation have separate licenses and attribution
under `samples/`. No new blanket license has been imposed on the supplied code.

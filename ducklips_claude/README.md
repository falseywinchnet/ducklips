# DuckLips speech engine

Status: **research prototype, working end to end** (2026-10-08). Not an
architecture decision and not an admitted GUI.Forms or File Manager component.

**GIVEN** (architect, 2026-10-08): do not synthesize speech; simulate a body
that is always speaking. Articulate the geometry and pliability of the
folds, throat and mouth; obtain its eigenstates so that the wave calculation
is a sparse solve that can be deformed anywhere; perturb it with tightening,
narrowing and noise; propagate the wave out of the mouth through air to a
microphone and sample that. A streamed behavioural cue sequence, built from a
dictionary and natural speaking behaviour, drives the body. Voices are
different bodies. Name: DuckLips.

DuckLips does that. Text becomes cues; cues move articulators; articulators
reshape a branched acoustic tube; the tube's eigenstates are re-solved by
warm-started sparse tree factorizations every millisecond; self-oscillating
vocal folds, airflow and turbulence drive those modes; the radiated flow
travels through air to two microphones. There is no recorded speech, no
formant filter bank and no learned model anywhere in the path.

See [DESIGN.md](DESIGN.md) for the mechanism and [FINDINGS.md](FINDINGS.md)
for measurements, failures and open problems.

## Build and run on the M4 Mini

From the Neo checkout (the Mini mirror is created by `m4build`):

```sh
/Users/ultimussecundai/.local/bin/m4build -- /bin/sh -c 'cd research/ducklips &&
  /opt/homebrew/bin/cmake -S . -B build -G Ninja -DCMAKE_MAKE_PROGRAM=/opt/homebrew/bin/ninja &&
  /opt/homebrew/bin/cmake --build build && ./build/ducklips_tests'
```

Speak:

```sh
./build/ducklips_say --voice alto --normalize \
    --text "Hello. My name is DuckLips." --out out/hello.wav
./build/ducklips_say --phonemes "HH AX0 L OW1 | AA1:1.2" --out out/phones.wav
./build/ducklips_say --voice baritone --room --text "In a small room." --out out/room.wav
```

Voices: `baritone`, `alto`, `child`, `gravel`, `whisper`. Useful flags:
`--report` (vowel fits and the phone timeline), `--probe FILE.csv` (per-chunk
physiology), `--trace START_S FILE.csv` (per-sample radiators),
`--diagnostics` (basis orthogonality every step), `--dictionary cmudict.txt`
(adds a CMU-format pronouncing dictionary), `--invariant
passive|charge|pressure`, `--room`.

Measure:

```sh
./build/ducklips_analyze out/hello.wav --spectrogram out/hello.png --clicks
./tools/measure_segments.sh baritone   # fricative and /h/ levels against /a/
```

GUI.Forms Audio integration (builds the pinned GUI.Forms audio service from
the `gui_forms` submodule; its CMake fetches and hash-checks miniaudio):

```sh
/opt/homebrew/bin/cmake -S . -B .build/gui-forms -G Ninja -DCMAKE_MAKE_PROGRAM=/opt/homebrew/bin/ninja \
    -DDUCKLIPS_WITH_GUI_FORMS_AUDIO=ON
/opt/homebrew/bin/cmake --build .build/gui-forms && ./.build/gui-forms/ducklips_gui_forms
```

`ducklips::GuiFormsVoice` is a `gui_forms::AudioGenerator`; the offline check
renders through `AudioEngine` and compares with a direct render. `--play`
opens the audio device.

Outputs go to `out/` and builds to `build/` (both ignored by Git).

## Layout

| Path | Contents |
|---|---|
| `include/ducklips/tree_operator.hpp` | tree-sparse symmetric matrix: shifted LDLᵀ, inertia count, solve |
| `include/ducklips/acoustic_network.hpp` | the airway as nodes (air volumes) and branches (air paths, walls, radiation) |
| `include/ducklips/modal.hpp` | eigenstate solver, modal resonator, energy-safe reprojection |
| `include/ducklips/sources.hpp` | two-mass vocal folds, Bernoulli orifice, intraoral pressure, noise |
| `include/ducklips/articulation.hpp` | articulators, the tract area model, constriction gestures |
| `include/ducklips/vowel_library.hpp` | formant targets inverted to tract postures through eigen-sensitivities |
| `include/ducklips/phonetics.hpp` | phone inventory, core lexicon, letter-to-sound fallback |
| `include/ducklips/planner.hpp`, `cues.hpp` | text to behavioural cue stream; lock-free cue queue |
| `include/ducklips/air_path.hpp` | radiators to microphones through air (and an optional room) |
| `include/ducklips/synthesizer.hpp` | the real-time engine (`render(std::span<float>)`, no allocation) |
| `include/ducklips/gui_forms_voice.hpp` | optional GUI.Forms `AudioGenerator` adapter |
| `tools/` | `ducklips_say`, `ducklips_analyze`, `ducklips_gui_forms`, calibration script |
| `tests/ducklips_tests.cpp` | numerical checks against dense references and closed forms |

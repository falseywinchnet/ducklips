# Fresh Papa render through the coupled engine

The current coupled source now has a full-word driver. It streams the prepared
mouth, lip and throat geometry together with lung pressure, effort and fold
rest-opening commands, preserving gas, tissue, radiation and passive-loss
state between intervals. No previous speech waveform is used.

The first attempt is `Papa.` with `bodies/coupled-oral-trial.body`: the reference
anatomy, a declared 1750 Pa lung-pressure setting, and no inhalation prelude.
The word is oral and all velum commands remain exactly zero. Source tension
is explicitly held at Q=1. Planned pitch is recorded but not applied; this
attempt does not complete pitch tracking, nasal speech or body calibration.
Those requirements remain part of the project.

## Control stream

`tongue_reach_audit --source-controls FILE.csv` exports the native scalar
trackers after applying planner cues and the body's rate limits. It includes
rest opening, lung pressure, intended pitch, effort, loudness and velum
opening. The optional export uses separate scalar state and leaves the
existing geometry, airway, lip and raw-cue outputs byte-identical. It never
exports a rendered source waveform or an imposed fold cycle.

The existing geometry exporter labels each block-end state with the block's
start time. The new CSV uses the same label for synchronization and records
the actual advanced time in `tracker_time`. The word driver linearly
interpolates the synchronized 1 ms endpoint values; the coupled controller
samples them at its actual trial midpoints. Lung pressure is multiplied by
the exported effort, clamped to the existing 0.3–2 range. Fold targets act
through the previously verified spring-bias forces.

`render_coupled_word.py` requires `--fixed-tension`, records the unapplied
pitch range, and rejects nonzero velum or varying playback-gain commands.
It validates frame synchronization and the geometry-preflight hashes, then
freezes its executable, body, source snapshot and inputs. Each successful
millisecond has a full continuation state, physical work/air checks and
input/output hashes. The live manifest records the active native PID.

The prepared input contains 1,366 frames. Native and independent geometry
checks pass all of them. The lip outlet closes at sampled time 0.532 s,
reopens at 0.670 s, closes at 0.742 s and reopens at 0.795 s. These are sampled
geometry checks; the acoustic solve still has to traverse every transition.

## Initial state and output size

The render starts from the prepared mouth pose at input time 0.399 s with
ambient gas, zero flow and folds at the commanded spring-rest opening. It
does not claim to have simulated the preceding 0.399 seconds. The first
pressure ramp begins in the following input interval. The final input frame
is 1.365 s, giving a 0.966-second requested render window.

`--initial-rest-equilibrium` initializes a fresh source at the spring-rest
command without the diagnostic's previous displacement perturbation. It
cannot be combined with restored physical state. In a zero-pressure held
mouth check, this source remains stationary; applying the following pressure
ramp passes whole-system work and air balance checks.

`--compact-memory` retains every intermediate site's energy and stored-air
observation, and saves its complete internal memory at the initial and final
states. It changes logging only. Full and compact outputs have exactly equal
physical observations and audit results. Four sanitizer runs cover the
equilibrium/ramp cases with both logging modes. Incomplete or contradictory
initialization flags are rejected.

## Active attempt and limits

The attempt is under `out/papa-coupled-oral-render`. Its `run.json` is the
authoritative live status. The committed `coupled_word_control_bridge.json`
binds input validation and a saved startup snapshot, not a completed word.
The current reference path remains far below the project's 10x-realtime
target. This attempt is for speech correctness; performance and perceptual
quality are not established by its numerical checks.

```sh
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/render_coupled_word.py out/papa-coupled-oral-inputs \
  --body bodies/coupled-oral-trial.body \
  --output out/papa-coupled-oral-render --fixed-tension
```

Do not restart that command while the recorded native process is live, or
overwrite its frozen files. Completion of the utterance, fresh audio export,
enhanced spectral inspection and human intelligibility assessment are still
pending. No diagnostic prefix is presented as a speech sample.

`actual_coupled_audio.py` now refuses a word-labelled manifest unless the
entire requested duration has completed. Running, failed and prematurely
labelled complete prefixes were all rejected before creating any audio
(`coupled_word_audio_guard.json`). Existing diagnostic manifests retain their
explicit diagnostic scope. A new word export preserves the intended word,
body hash, applied controls and unapplied pitch in its receipt.

The following completion process is waiting for this attempt. On successful
completion it checks every trajectory hash and physical join, exports the
new 48 kHz waveform with one fixed gain selected from its whole peak, and
renders the full enhanced field in original color at native pixel resolution.
It checks the WAV against the exported physical signal and the PNG against
every enhanced-field pixel. It emits no sample if the word render fails.
Positive full-word export and listening remain pending.

```sh
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/finish_coupled_word.py out/papa-coupled-oral-render/run.json \
  --out out/papa-coupled-oral-listening --wait
```

A read-only audit of the first 122 completed segments verified all 121 joins:
maximum source-energy reconstruction difference 1.90e-19 J, maximum whole
work residual 1.80e-15 J and maximum air residual 3.84e-20 m³. This validates
the completed preparation prefix, not the remaining speech or its quality.

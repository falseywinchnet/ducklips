# Continuous voice trajectory to radiation audio

`tools/coupled_phonation.py` connects completed native coupled-source
trajectories to a longer physical continuation and an audio export. This is
the eight-section fixed-tube control at 1400 Pa lung pressure. It is not the
articulated NATO renderer, and a successful export is not evidence of realistic
speech or settled voicing.

Every continuation restores the fourteen source values, eight tract charges,
eight tract momenta, radiation momentum, and 192 thermal plus 192 viscous
memories. The driver checks the restored source, gas mass and energy; the
tract and all loss memories must match exactly. Adjacent traces share one
endpoint, which is retained once. No waveform cycles are repeated. Native
failures retain their input, partial output and stderr and stop the driver.

The executable is checked against its build receipt before work and before
each segment. A new output directory prevents overwriting a live run. Each
completed segment records its command, input/output hashes, wall time and
restart discrepancy. Five-millisecond checkpoints allow further continuation
from the latest completed state. They do not avoid the computational cost of
the current accuracy-reference engine; production efficiency is unfinished.

## What the audio represents

The existing aperture boundary defines pressure and radiated power by

```
P = R (Q_outlet - Q_radiation)
power = P² / R
```

The exported signed signal is `sqrt(R) * (Q_outlet - Q_radiation)`, in
`sqrt(W)`. Its square is that boundary's instantaneous radiated power. This
uses the actual tract outlet and radiation memory, not a fold waveform or
interior pressure substituted for emitted sound. It remains a radiation
power signal: microphone pressure, distance, directivity and calibrated SPL
are not assigned to it.

Uniform native observations are resampled to 48 kHz with a polyphase
anti-alias filter (Kaiser beta 8.6). The final endpoint provides filter context;
playback contains only whole output intervals within the simulated duration.
The exporter saves the original observations and resampled values as float
arrays. The playable WAV is PCM16 with a declared fixed gain, no peak
normalization, no repetition and no time stretching. A clipping gain is
rejected. Exterior filter padding is zero; endpoints of a short excerpt can
therefore include filter transients. The recorded trapezoidal power integral
is an observation integral, not the native integrated energy receipt.

## Checked behavior

The [control receipt](coupled_phonation_controls.json) records:

- A signed analytic 1 kHz radiation signal has the expected energy and phase.
  Resampling changes its interior amplitude by at most 1.56e-7 in this fixture.
- A 35 kHz component sampled at 100 kHz is attenuated before the 48 kHz output;
  its residual RMS is 1.33e-5 of its original RMS in the checked interior.
- A real 20 µs continuation of the completed 5 ms onset preserves gas mass
  exactly, source state within 2.17e-15 in the declared scaling, and all tract
  and loss memories exactly. The combined 503 native samples have one shared
  boundary and produce 240 complete playback intervals.
- Its outlet and radiation observations match the original long continuation
  exactly at those three common sample times. Source pressure differs by at
  most 2.17e-7 Pa and tract pressure by 6.41e-10 Pa. Segment-local time
  arithmetic differs by a few floating-point ulps.
- Irregular observation clocks, negative radiation resistance, reset thermal
  memory, a changed source effort and incomplete native segments are rejected.

The first driver test exposed a decimal timing issue: subtracting 5 ms from
5.02 ms produced a value just above 20 µs, and the native `ceil` sampling rule
created three observations instead of two. The driver now decimalizes the
segment duration and verifies the resulting observation spacing before
launching. It also rejects irregular concatenated times before exporting.
The failed experiment is retained under `out/coupled-phonation-driver-control`;
the successful check uses the separate `-v2` directory.

## Use

Use the scientific Python environment. Supply completed contiguous seed
segments from the same native body, drive and executable; source-state checks
cannot infer an unrecorded change in body parameters.

```sh
python tools/coupled_phonation.py run \
  --binary out/joint-contact-v2/probe \
  --build out/joint-contact-v2/build.json \
  --seed out/joint-contact-v2/phonation-start.jsonl \
  --seed out/joint-contact-v2/phonation-continued.jsonl \
  --out out/coupled-phonation-long --duration .25
```

The second seed is usable only after its native run has completed. The
driver exports the entire accumulated radiation trajectory after every
successful segment. To export an existing completed set independently:

```sh
python tools/coupled_phonation.py export \
  out/joint-contact-v2/phonation-start.jsonl \
  --out out/coupled-onset-radiation
python tools/coupled_phonation_controls.py
```

A few milliseconds of exported onset are a tooling check, not a listening
sample. Longer pressure-driven motion remains in progress. No new NATO
pronunciation claim follows from this work.

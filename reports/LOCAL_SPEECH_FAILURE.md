# Rejected speech and inlet correction, 2026-10-10

The user rejected every recording from `c528814a4877caca441518e68fd5b771642dd8bf`.
Papa was high-pitched noise, a cut-off first sound, then “pa”; November was
“click-click-ember”; Foxtrot was only noise; the alphabet was clicks and garbage.
The first two associations follow presentation order. Exact feedback and played
WAV hashes are in `engine_recovery.json`. None is an approved pronunciation.

This investigation does not recover intelligible speech. It corrects one source
mechanism and separates it from unresolved geometry and onset failures. No new
alphabet was rendered, and diagnostic WAVs were not offered as listening results.

## Isolation

The rejected executable was retained. The original Papa trajectory was traced
with unchanged body, dictionary and controls. Its first vowel occupies
0.678–0.731 seconds. Only two glottal-area peaks with prominence at least
2e-6 m² occur there, at 0.716063 and 0.724417 seconds. Fold-area peaks alone do
not establish an audible vowel; the next /p/ is already approaching closure.

Disabling both turbulence sources preserves the early high-frequency output.
For a stronger check, a diagnostic copy of the body disables audible breathing,
and CLI overrides set exhalation pressure and both turbulence gains to zero.
This removes inhalation as a possible explanation too. The moving geometry
still produces substantial output:

| Papa interval | Rejected engine lip acceleration RMS | Geometry-only RMS |
|---|---:|---:|
| 0.400–0.600 s, before first phone | 0.018871 | 0.020109 |
| 0.678–0.731 s, first vowel | 0.052908 | 0.030525 |
| 0.850–1.050 s, later voiced region | 0.074232 | 0.003432 |

Units are m³/s², before the microphone path or playback gain. These are source
diagnostics, not perceptual scores. Geometry motion can legitimately radiate;
this test identifies the source of the unwanted output, not a requirement that
all articulatory motion be silent. Reducing the original control block from
48 to 12 samples lowers early RMS to 0.006228, demonstrating strong sensitivity
to its temporal representation. That setting was not adopted as a fix.

The enhanced Papa field completed independently. Its source hash matches the
rejected playback exactly. `papa-enhanced-native-portrait.png` rotates the full
684 × 1371 field without interpolation or cropping and retains the original
magma mapping, reference and dynamic range. Time increases downward and
frequency increases rightward. No ordinary STFT comparison was used.

## Source correction

The old inlet was a memoryless Bernoulli resistance integrated separately from
the tract. Charging its small inlet compliance with a midpoint step can cross
the reservoir pressure. The pressure sampled by the folds then depends strongly
on the split-update phase. An experimental exact resistive charging update
equilibrated that isolated inlet each step and stopped phonation entirely; it
was rejected rather than retained as a supposed stability improvement.

The inlet now carries the momentum of the air crossing the folds. With opening
area A and fold thickness L, its inertance is M = ρL/A. L comes from the body's
lower plus upper fold thickness. The state obeys

```
dp/dt = P_lung - P_inlet - β Q |Q|,  Q = p/M,  β = ρ/(2 A²)
dq_inlet/dt = Q
```

One local midpoint exchange has a scalar quadratic solution, evaluated in a
cancellation-free form. This adds one persistent momentum, no global solve,
iteration, adaptive step, prescribed oscillator or audio filter. Source kinetic
energy is included in the reported energy. At an exact/numerical seal the
passage momentum is cleared and retained cavity air is unchanged.

For Papa's first-vowel interval, inlet pressure changes from −55…870 Pa to
−144…169 Pa. Three fold-area peaks occur at 0.710229, 0.718688 and 0.727 seconds.
This supports the coupling correction; it does not prove the first syllable
works. Early geometry-driven noise persists and is slightly larger after the
source change. No quality improvement is claimed from this render.

## Checks and remaining work

Release and address/undefined-sanitizer local tests pass. New isolated tests
check reservoir-work passivity, supplied-air conservation, reversed breath,
kinetic energy accounting, closure and timestep refinement at the same elapsed
time. Existing local streaming checks still establish identical chunked output,
no render allocations, and no nonfinite-sample substitution or clipping in their
test trajectory. They do not validate moving-source work or intelligibility.

The next necessary representation repair is in moving geometry: the current
upper-gap midpoint sections and compliance floor are being treated as physical
air volumes. They omit lower-tongue air and can change abruptly when the blade
crosses a section. Piecewise-linear geometry also introduces velocity changes
at control boundaries. A correct repair must carry continuous occupied volumes
and their swept air, including contact, rather than filter the emitted clicks
or switch off displacement pressure. The onset schedule must then be checked
against the actual fold and lip trajectories so that an intended vowel has time
to exist. The present correction leaves both tasks unresolved.

Diagnostic data, exact commands, hashes and the new executable/source snapshot
are retained locally in `out/local-failure-investigation/run.json`.

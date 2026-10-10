# Rest-opening actuation in the current coupled engine

The actual-mouth adapter now accepts `--fold-rest-target START END` in square
metres. The values are signed spring-rest targets: a negative value commands
compression. The input's geometrical rest area stays fixed. Forces move the
existing tissue, and the shared solve evolves gas, momentum and passive memory
without resetting them. This does not prescribe a fold cycle or render speech.

For baseline area A0, target A, fold length L, layer stiffness k and the already
configured tension/asymmetry factor Q, each side receives

```
displacement target d = (A - A0) / (2 L)
external force F = k Q d
```

This is a spring-bias actuation model, not a calibrated laryngeal muscle model.
Both layers receive the same displacement target, so their coupling spring
needs no additional bias. The core already includes external force in the
implicit equations and records its work as the sum of force times tissue
displacement. Commands are sampled at each adaptive trial's true midpoint,
including half steps and retries. Merely changing the command leaves the
geometrical baseline and current state intact. The command has no additional
smoothing; the caller must supply the intended target trajectory.

The Python actual-mouth runner exposes the same optional pair through
`fold_rest_target`. Existing callers retain their old behavior. A fresh
unflagged replay of the previously saved short Papa vowel load matches its
entire parsed output exactly.

## Verification

`tools/fold_rest_control_reference.py` uses the freshly prepared Papa vowel
geometry at 0.900 seconds, holds the mouth fixed, and sets lung pressure to
zero to isolate actuator work. It runs 40 microseconds at a constant 4e-6 m2
rest target, then continues the complete state for another 40 microseconds
while returning the target to the 3e-6 m2 baseline. All six runs at local
tolerances 1e-4, 2.5e-5 and 6.25e-6 complete.

The constant-force work has a separate exact endpoint integral F dot (x1-x0).
Its discrepancy from recorded external work is below 2e-24 J. The whole-system
energy residual is below 1.1e-20 J and air balance residual below 2e-20 m3.
Successive final pressure differences decrease from 0.00375 to 0.00186 Pa;
final tissue-position differences decrease from 5.46e-12 to 4.40e-12 m.
These are refinement comparisons, not errors against an independent complete
trajectory solution.

The final signed-target implementation reproduces the positive-target run
exactly. A -1e-6 m2 command also passes the work and air checks. Sanitizer builds
pass both constant-force and continued-ramp cases, including the work audit.
Malformed/nonfinite commands and use without the coupled accuracy controller
are rejected. The frozen run, source snapshot, native inputs and final-build
verification are under `out/fold-rest-control`; the compact committed receipt
is `fold_rest_control.json`.

## Remaining speech work

This connects one missing physical control. It does not yet consume the full
planner cue stream, apply pitch/tension commands with their associated work,
or add the missing nasal branch. Sustained source closure also remains under
investigation. There is no newly rendered word or new listening score from
this change.

The separate near-closure comparison in `source_update_comparison.json` tests
an optional dissipative update from the same saved state. Over 0.1 microseconds,
the fitted update uses 128 accepted leaves; the dissipative update uses 2,048,
and 4,096 at the tighter tolerance. Numerical loss approximately halves with
refinement, but this does not establish a closure fix. The experimental API
was removed; its executable, source snapshot and outputs remain preserved.

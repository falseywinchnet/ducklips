# Lung-pressure control without a source reset

The actual-mouth adapter accepts an optional `--reservoir-end PA`. The input
reservoir pressure is the start of a linear control over that geometry
interval. `planar_continuous_reference.run(..., reservoir_end=...)` exposes
the same control to the stateful Python driver. The option requires the
continuous source and coupled accuracy controller.

Each trial samples pressure at its own physical midpoint. Full trials,
rejected trials and the two accepted half steps therefore sample the same
time-dependent control consistently. This is midpoint integration of a
changing input, not an exact integration claim. The existing native gas law
converts pressure to reservoir mass effort, and the work receipt accumulates
that effort against transferred reference volume on accepted steps. Source,
graph, radiation and loss-memory state continue across chunks.

The [validation](reservoir_control_validation.json) starts from the completed
50 ms actual-mouth state. A deliberately accelerated 1400-to-1750 Pa ramp
lasts 0.2 ms in the held final mouth. It checks:

- The default constant-pressure output matches the preceding frozen native
  executable exactly. An explicit equal-pressure endpoint also gives exactly
  the same physical output.
- Splitting the ramp at 0.1 ms retains the state and reproduces its endpoint
  to numerical precision.
- Independent gas-law bounds enclose each observation interval's reservoir
  work. Whole-trajectory energy and air balance remain checked.
- The successive endpoint pressure differences are 0.00817, 0.01523 and
  0.000657 Pa as maximum trial size decreases from 10 to 2.5 to 0.625 to
  0.3125 microseconds. The early adaptive patterns are not monotone. The last
  two runs have no rejected trials and their difference decreases substantially.
  This is local refinement evidence, not a global accuracy guarantee.
- An [ASAN run](reservoir_control_asan.json) completes with no sanitizer
  diagnostic; its largest sampled pressure difference from Release is
  6.79e-10 Pa.

The ramp is a numerical control fixture, not a physiological breath or a
new speech recording. It supplies a physical input needed by streamed speech;
it does not prescribe the fold oscillation or the output waveform.

## Sustained source remains unresolved

The separately frozen fixed-tube continuation has terminated at its trial
limit. It completed 25 ms with a restartable endpoint, then recorded another
3.95 ms before failing within the next observation interval, at approximately
28.954497 ms total. The [failure receipt](coupled_phonation_long_failure.json)
preserves its exact executable identity, trace, input and diagnostic.

At the last emitted observation the lower source opening is closed, the upper
opening is about 1.17e-7 m2, and the first stored gas mass is 6.67e-14 m3 of
reference volume. Its corresponding effort is about 606 kPa. The controller
reports `JointAccuracyFailure::trial_limit`; this observation does not establish
whether additional trials would resolve the compressed gas/contact state.
Neither the tolerance nor the trial budget has been changed to hide the failure.

The independent actual-mouth 1750 Pa run remains separate. Its first 22 ms
have completed as of this record; its frozen executable is unaffected by
the new optional pressure control.

## Next audio deliverable

The user clarified on 10 October that audio deliverables must be newly
rendered speech from the engine under development. Previously supplied
spatial-pressure recordings and short onset diagnostics do not meet this
request. No new word or intelligibility improvement is claimed by this work.

Reproduce the pressure-control validation with:

```sh
cmake --build build --target ducklips_planar_inertive_motion_probe
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python tools/reservoir_control_validation.py
```

The recorded seed and articulation fixtures are required. Native intermediate
files are preserved under `out/reservoir-control/`.

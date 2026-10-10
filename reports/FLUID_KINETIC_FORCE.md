# Differentiate the fluid inertia before adding constant tissue mass

The moving-wall force calculation was losing small fluid-inertia changes
inside a matrix dominated by constant tissue masses. The source previously
differenced inverse total-inertia matrices. A first comparison that
differenced the total inertia instead still included the large constant
diagonal terms and retained much of the numerical noise.

For kinetic energy K = pᵀG⁻¹p/2, the derivative at fixed canonical momentum is

```
v = G⁻¹ p
dK = -vᵀ(dG)v / 2
```

Only the fluid part of G depends on the two gas moments and two opening
coordinates being differentiated. The tissue masses and attached outlet
column are constant in this derivative. The source now differences the fluid
inertia before adding those constants, maps that derivative into the coupled
coordinates, and contracts it with the midpoint velocities corresponding to
the existing endpoint momenta. The same identity applies to the dissipative
update. The discrete work correction, physical loss laws, contact rule and
coupled error tolerances remain in place.

## Independent moving-wall check

The earlier small-gas force fixture held the walls still. This check uses the
saved 28.954 ms fixed-tube state, including both folds' actual velocities.
The independent calculation integrates scalar axial kinetic energy at fixed
physical velocity and negates its coordinate derivative. It reconstructs
the pressure field with separate gas quadrature and solves for fixed gas
moments when varying the geometry. It does not invert the native inertia
matrix or add tissue constants before differencing.

At the default difference step, the first gas-moment derivative has these
errors against the independent result:

| Calculation | Absolute error, Pa |
| --- | ---: |
| Inverse total inertia | 0.01749 |
| Direct total inertia | 0.003356 |
| Direct fluid inertia | 8.31e-7 |

Across three matched difference sizes, the fluid-only errors remain below
8.31e-7 Pa for gas-moment derivatives and 1.05e-13 Pa m for opening-area
derivatives. The independent gas derivative changes by less than 2e-6 Pa
between its final two refinements. These are checks of this model's force,
not physiological calibration or a complete trajectory oracle.

The raw coupled full-step/half-step comparisons at 10, 25, 100 and 500 ns now
all solve from that saved state. Previously the 25 ns half-step and the
100/500 ns full steps failed. Solving a trial does not establish its accuracy;
the coupled controller still has to accept its measured local error.

Six continuous-source/coupled/contact regressions pass in Release and in the
sanitizer build. `tools/moving_source_force_reference.py` reproduces the
independent force calculation, while `continuous_force_probe` retains the
inverse and total-metric alternatives for comparison.

## Closure remains incomplete

A further run starts at the fully saved 28.956 ms state. The updated source
captures another 0.740 microseconds with a positive upper opening of
6.9509e-10 m2, then fails at the minimum permitted time step. The two gas
moments are 2.0107e-20 and 2.6308e-16 m3. The updated force improves numerical
precision but does not resolve disappearance of this remaining gas pocket.
The completed prefix retains checked work and mass balance; its maximum
residuals are 5.282e-16 J and 1.336e-21 m3.

`joint_contact_probe --checkpoint-every N` now optionally saves the complete
physical state at intermediate observations. A comparison with and without
that flag produces identical physical observations. This preserves the last
successful state when a subsequent interval fails, avoiding a replay solely
to recover it. Default checkpoint behavior is unchanged.

An earlier attempt to bound the Newton perturbation by each remaining gas
moment did not improve the initial trial-failure pattern and was not adopted.
Its frozen source and instrumented failure traces remain under
`out/source-local-jacobian`. They locate the failure in Newton's line search.
The total-metric and fluid-only experiments are separately preserved under
`out/source-metric-gradient` and `out/source-fluid-gradient`.

## Completed actual-mouth source run

The independently running 1750 Pa actual-mouth diagnostic completed 80 ms
using its executable and source snapshot frozen at launch. It has 79 verified
state joins, positive source gas throughout, and zero source closures or
reopenings. Oscillation continues through the final window. Its largest
per-segment work and mass residuals are 7.452e-16 J and 3.153e-20 m3.

This trajectory includes one millisecond of mouth movement and otherwise
holds the mouth. Its amplitude and peak intervals vary; it does not prove a
stable periodic voice. It is a usable saved state for further control work,
not a word recording or a regression of the newly changed force code. The
completed receipt is `actual_coupled_1750_completed.json`. No audio from this
diagnostic is offered as newly rendered speech.

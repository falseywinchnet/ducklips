# Coupled motion of continuous gas and tissue

`ContinuousGlottis` now advances the continuous gas field and four tissue
masses together. It uses the shared fluid/tissue inertia, gas compression,
elastic contact, tissue damping and finite BGK slit loss. Pressure drives the
motion; no source waveform or tissue oscillation is prescribed.

This is a native mechanical reference within a fixed port topology. It can
advance an open source or an already sealed lower/upper cap with remaining
wedge-shaped air. It does not yet carry a trajectory through sealing or
reopening, vary the rest opening, or connect this source to the speech tract.
Those are required next steps. No new speech recording or pronunciation score
comes from these controls.

## One kinetic energy for the moving air and tissue

Let J0, J1, J2 denote accumulated inlet, internal and outlet reference-volume
transfer, and let x contain the four tissue positions. The coordinates are

```
q = [J0, J1, J2, x_lower_left, x_lower_right, x_upper_left, x_upper_right]
M0 = initial_M0 + J0 - J1
M1 = initial_M1 + J1 - J2
a = rest_area + width*(x_lower_left + x_lower_right)
b = rest_area + width*(x_upper_left + x_upper_right).
```

J1 corresponds to the continuous field's integrated internal flux, not a
flow through the old chamber seam. The linear rate transformation T gives

```
w = T qdot = [J0dot, J1dot, J2dot, adot, bdot]
G = T^T L T + tissue_mass_matrix
kinetic_energy = qdot^T G qdot / 2
p = G qdot
Hamiltonian H(q,p) = p^T G^-1 p / 2 + gas_energy + tissue_springs + contact.
```

The gas state is recovered from its two overlapping mass moments at the
current contour. G therefore depends on both mass and geometry. Its
configuration derivatives produce real momentum/force terms; retaining only
pressure forces with separate passage inertias would omit them. Fluid motion
induced by the walls remains part of G, including its cross terms.

A closed exterior port is removed from the admissible coordinates. Its flow
and reported momentum are zero; the remaining coupled metric is inverted.
This is phase-wise constraint elimination, not a rule for discarding momentum
when an open cap closes. The latter transition is not implemented. Trial steps
that change port topology are rejected without changing the accepted state.

## Work-consistent time update

The implementation uses dimensionless coordinates internally for matrix
conditioning; external values remain SI. It solves for the seven coordinate
increments, eliminates the final momentum using the symmetric endpoint
inverse metrics, and solves the resulting nonlinear force balance.

For the active coordinates, define

```
qdot_bar = (G0^-1 + G1^-1)(p0 + p1)/4
q1 - q0 = dt*qdot_bar.
```

The configuration part of the finite kinetic change is

```
[ p0^T (G1^-1-G0^-1) p0 + p1^T (G1^-1-G0^-1) p1 ] / 4.
```

A symmetric midpoint derivative plus a discrete-gradient correction accounts
for this change and gas/contact energy exactly to solve/rounding tolerance.
The correction lives in the four actual configuration variables `[M0,M1,a,b]`.
It therefore preserves the cyclic common-volume-transfer coordinate. Tissue
springs use their exact midpoint gradients. Kinetic derivatives use centered
configuration differences; their accuracy is checked independently below.
This remains a numerical mechanical reference, not an optimized renderer.

At the midpoint configuration, the positive slit-loss matrix and tissue
damping act on the same average velocity. Reservoir pressure is converted to
its gas enthalpy effort. The outlet accepts a free effort plus a nonnegative
linear load slope. An explicit optional outlet kinetic-energy loss coefficient
(default 1) gives the additional effort

```
coefficient * rho0 * Qout*abs(Qout) / (2*r_out^2*A_out^2).
```

That is a stated exit-loss closure, not a measured anatomy-specific result.
The scalar work ledger is

```
Delta H + fluid_loss + tissue_loss + outlet_work
        - reservoir_work - external_tissue_work = residual.
```

No numerical damping term or energy deletion is inserted. Rest opening is
held fixed in the present advance API, so there is no unaccounted actuator
work. A future changing-rest implementation must include that work explicitly.

## Precision needed near ambient pressure

The first equilibrium/derivative test exposed a failure in inverse mass
recovery near ambient conditions. The gas calculation converted a tiny density
departure through `1 + departure`, then subtracted one again to recover the
energy. That lost enough relative information to defeat the storage quadrature.

`GasStorage::evaluate_effort` now computes logarithmic density, pressure and
compression energy directly from enthalpy effort using `log1p`/`expm1` and the
existing small-departure energy series. The density-based API retains its
behavior at very small positive density. No pressure or density floor is
introduced. [Fresh 80-digit checks](continuous_glottis_storage_audit.json)
cover the original 72 configurations plus six signed effort fields as small
as 1e-12 Pa. The largest scaled storage discrepancy is 1.43e-15; the largest
relative energy discrepancy in those six near-ambient cases is 4.34e-16.
The previously failing inverse-moment case now succeeds.

## Accuracy is checked separately from conservation

[The motion audit](continuous_glottis_reference.json) uses four nonlinear
controls: inviscid free motion, 700 Pa reservoir drive with loss, sealed lower
cap, and sealed upper cap. Each includes independently moving left/right
masses with asymmetric tissue parameters. The sealed controls retain positive
air storage while their closed exterior flow remains zero.

The coarse 100/200/400-step trajectories over 100 microseconds conserved
energy and mass but failed the accuracy criterion. [Those failures are
retained](continuous_glottis_initial_accuracy.json), including the additional
lower-cap refinement failure at 3200 steps. The grid was refined rather than
the acceptance threshold relaxed. The final control grids and measured
convergence ratios are retained in the motion audit.

The accepted refinement comparisons are:

| Control | Steps over 100 microseconds | Endpoint difference ratio |
| --- | --- | --- |
| Free, inviscid | 800 / 1600 / 3200 | 3.609 |
| Driven, dissipative | 800 / 1600 / 3200 | 4.046 |
| Lower cap sealed | 1600 / 3200 / 6400 | 4.660 |
| Upper cap sealed | 800 / 1600 / 3200 | 3.790 |

These ratios approach four for a second-order method and pass the original
3.3–4.7 interval. They measure differences in normalized endpoint mass,
position, flow and velocity; they do not bound every trajectory value or
establish accuracy of the longer exploratory run.

There is also an independent linear limit. Separate continuous-field
quadrature constructs the ambient fluid metric. The gas stiffness follows
analytically from compliance C and geometry derivative D:

```
du = C^-1 (dM - D da)
K_gas = [B, -D W]^T C^-1 [B, -D W],
```

where B maps the three volume transfers to mass moments and W maps tissue
positions to the two control areas. With tissue springs and masses included,
a matrix exponential solves the resulting 14-state linear ODE directly.
The native nonlinear motion approaches this solution as perturbation amplitude
and step size decrease. At amplitude 1e-4 and 1600 steps over 10 microseconds,
its maximum normalized endpoint error per unit amplitude is 3.95e-6, passing
the original 1e-5 criterion. This normalization scales masses by 1e-8 m³,
positions by 1e-4 m, flows by 1e-5 m³/s, and tissue speeds by 0.1 m/s; it is
not a blanket relative error on every field.

The linear calculation identifies tissue modes near 108, 127, 180 and 211 Hz,
and short-column pressure modes near 63.6 and 140.3 kHz. These fast modes
explain why energy conservation did not establish trajectory accuracy on the
coarse grids. They also identify a physical timescale that a future efficient
representation must address. Performance optimization remains parked.

A separate [nonlinear force audit](continuous_glottis_tangent_reference.json)
uses substantial airflow and moving tissue. Independent quadrature and mass
recovery differentiate the fluid metric at fixed mass/area. The expected
canonical force includes `-v^T (dG) v / 2` alongside gas and tissue forces.
Three decreasing native steps approach that independently computed initial
momentum derivative. The refined reference derivative changes by 2.78e-10 in its stated force
scaling, and the native initial-rate Richardson estimate differs by 3.02e-8.
This tests the nonlinear inertia terms that an ambient linear-mode test would
not exercise.

An additional 1 ms, 700 Pa run advances without loss of mass or work balance
and develops pressure-driven tissue motion. It is an exploratory longer run,
not a time-refined accuracy result, a completed contact cycle, or evidence of
self-sustained voicing.

[Validation receipts](continuous_glottis_validation.json) record five focused
CTests and address/undefined-behavior sanitizer tests, with identical release
and sanitized outputs for four short motion controls and the nonlinear-force
probe. The earlier independent loss fixture is reused with its input/output
hashes verified; inverse-moment and prescribed-geometry seam checks are fresh.
The complete source/binary/output hashes are recorded in the audits.

## Reproduce and continue

Build the coupled probe and existing constitutive probes, then run:

```sh
OPENBLAS_NUM_THREADS=1 python tools/continuous_glottis_experiment.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_glottis_reference.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_glottis_tangent_reference.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_glottis_storage_audit.py
# Requires corresponding .build/asan targets configured with
# -fsanitize=address,undefined -fno-omit-frame-pointer:
OPENBLAS_NUM_THREADS=1 python tools/continuous_glottis_validation.py
```

The next physical/numerical task is carrying the coupled momentum through
cap sealing and reopening without erasing air, momentum or energy. Sustained
phonation, connection to the moving tract, streamed body/control changes and
new NATO listening remain unverified. The local BGK loss, simple medial
contour and lumped tissue/contact assumptions also remain model limitations.

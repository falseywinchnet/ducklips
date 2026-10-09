# Pressure-driven onset in the continuous source

The full-span continuous gas/tissue model has a pressure-driven oscillatory
instability in its symmetric tissue subspace. A separate implementation of the
continuous equations predicts decay at 700 Pa and growth at 1400 Pa. Native
small-perturbation runs approach both predictions under time refinement. This
establishes an onset mechanism in this model; it does not establish a stable
large-amplitude voice, physiological calibration, or intelligible speech.

The investigation also found and fixed a solver failure near steady throughflow.
The fix resolves a numerical force error without changing the configured tissue,
gas, pressure, geometry, or constitutive loss laws.

## Cancellation in the configuration work correction

The existing discrete gradient starts with a midpoint force, then adds a
correction parallel to the actual mass/area change `dz`:

```
correction = (endpoint energy change - midpoint force dot dz) / (dz dot dz)
gradient += correction * dz.
```

Near a steady flowing state, both endpoint energies remain finite while `dz`
can be extremely small. The difference of endpoint energies contains absolute
floating-point error. Dividing that error by the displacement creates a noisy
force. The tiny-perturbation control exposed this: a normalized state increment
with squared length about `1.44e-15` produced energy discrepancies of order
`1e-15` to `1e-14` in the normalized energy. The resulting force correction
prevented Newton convergence despite physically small, smooth motion.

The solver now applies the endpoint correction only if its discrepancy exceeds
an endpoint-arithmetic scale:

```
64 * double_epsilon * (
    abs(U0) + abs(U1)
    + componentwise absolute kinetic quadratic-form scale
).
```

The kinetic scale includes both endpoint inverse metrics and the same momenta
used by the selected update. Below that scale the midpoint gradient is retained.
The actual energy residual is still reported. `energy_roundoff_scale` and
`roundoff_limited_gradient` expose when this occurs. The scale is an arithmetic
error estimate, not a certified bound on all quadrature/inverse-solve errors.
Newton's force tolerance has not been increased.

This is relevant to the dynamics: machine roundoff should not become a physical
force just because a body is close to equilibrium. The 700 and 1400 Pa CTest
controls exercise the previously failing regime. The longer low-amplitude
controls check its dynamics against an independent matrix exponential, rather
than merely requiring that it no longer stops.

## The cyclic transferred-volume coordinate

A second failure appeared only after long throughflow. Reconstructing the same
physical state with `prepare` made a failed step succeed. The native coordinates
were storing cumulative inlet, middle and outlet transfers, then subtracting
these large values to recover small gas moments. Their common transferred-volume
coordinate is cyclic: the configuration depends only on the two differences.

After an accepted step the solver now retains the two target mass moments as
its new origins and resets all three transfer positions to zero. Canonical
momenta, tissue positions, reservoir work and the step's transported volumes
are retained. This is a coordinate translation, not a reset of the physical gas
or its momentum. It avoids losing pressure/Jacobian precision as an irrelevant
common transfer grows. A continuous-versus-restarted 10 ms control checks the
physical invariance. The long pressure-driven run now passes its earlier stop
and enters contact and reopening.

## Separate continuous equations

`continuous_phonation_reference.py` independently implements the open full-span
profile with Gauss integration, local continuity, the checked BGK conductance,
and symmetric paired tissue masses. It does not call the native time update.
Its physical state is

```
y = [mu0, mu1, a, b, Qin, Qmiddle, Qout, adot, bdot].
```

The two mass-moment derivatives give

```
mudot = C^-1 ([Qin-Qmiddle, Qmiddle-Qout] - D [adot,bdot]).
```

The shared metric is the separately integrated fluid metric plus effective
paired-tissue masses `m/(2 width²)` in the area coordinates. With `v` containing
the three fluxes and two area rates, the canonical-to-velocity conversion is

```
p = G v
G vdot = canonical force - Gdot v.
```

Metric derivatives are evaluated with a complex perturbation of the smooth
open-domain formulas. They are transformed from enthalpy/area coordinates to
mass/area coordinates using `C` and `D`, retaining the same kinetic shape force
as the Hamiltonian formulation. The independent matrices agree with native
quadrature. A substantial-flow, moving-wall control compares its instantaneous
rates with three decreasing native time steps and a Richardson extrapolation.

Stationary throughflow has `Qin=Qmiddle=Qout=Q` and zero area rates. The five
stationary momentum equations are solved for `mu0,mu1,a,b,Q`. Linearizing the
nine-state ODE then supplies growth rates and frequencies. The onset criterion
is a complex mode crossing from negative to positive real part. This use of a
linearized pressure equilibrium follows the stability framework discussed in
[Zhang (2011), *On the difference between negative damping and eigenmode
synchronization as two phonation onset mechanisms*](https://www.surgery.medsch.ucla.edu/spl/papers/2011JASA04.pdf).
That paper does not validate this engine's geometry or predicted threshold.

## Checked onset, with numerical damping separated

For the configured rest area `3e-6 m²`, default tissue parameters and outlet
loss coefficient one, the independent model gives:

| Reservoir pressure | Dominant growth rate | Frequency |
| --- | ---: | ---: |
| 700 Pa | −10.5734 /s | 146.4517 Hz |
| 1400 Pa | +19.7552 /s | 148.7651 Hz |

A pressure continuation places the model's onset near **870.767 Pa**, at
**147.3203 Hz**. Changing Gauss order from 24 to 48 and halving the linearization
increment shifts this pressure by about `0.00022 Pa`. This precision describes
numerical agreement within the model, not uncertainty in a human voice.

Native runs use equal-and-opposite small modal perturbations to cancel common
equilibrium drift. Their difference is compared with the independent matrix
exponential after two milliseconds. Three step sizes, 10, 5 and 2.5 microseconds,
reduce endpoint error by approximately a factor of two each time. Extrapolated
native growth rates agree with the continuous reference. The dissipative update
adds appreciable damping at finite time steps; this cannot be counted as tissue
or aerodynamic loss. Full numbers and hashes are in
[the independent onset receipt](continuous_phonation_onset.json).

## Large motion is a separate requirement

The nonlinear investigation uses steady 1400 Pa reservoir pressure, zero
specified downstream effort and no external tissue force. Earlier startup runs
showed growing motion, then stopped near 118 ms while the passage was still
open. Their failed outputs are retained. A continuation from a saved state,
using the roundoff fix, passes that earlier stopping region.

However, the dissipative update's coarse large-motion run appears to settle
while a finer run keeps growing. The coarse appearance alone is therefore not
evidence of a converged limit cycle. An independent stiff ODE integration of
the open-domain equations supplies an additional trajectory reference. It
stops when either area first reaches `1e-8 m²`; this is an explicit diagnostic
stop before contact, with no area floor or modification of the equations.
The native symmetric update is also compared while the passage remains open,
so continuous dynamics can be distinguished from damping introduced by the
all-step dissipative update.

The experimental driver now uses the symmetric update within a contact phase,
and accepts a dissipative step only when an exterior port actually changes
between open and closed. If neither trial succeeds, the driver bisects the
step, with a bounded recursion depth. Trial copies and per-substep work ledgers
make rejected attempts transactional. This selects a numerical update from the
computed geometry; it does not schedule closure or prescribe an oscillation.

From a common pre-contact state on the independent open trajectory, steady
1400 Pa pressure produces three closures and three reopenings in the subsequent
24 ms. The contact experiment refines the maximum step, records every accepted
substep in the work/mass totals, and verifies exactly zero exterior flow in
sampled closed states. The cycles are transient, and the endpoint differences
near the next closure still matter; a decreasing difference does not establish
a perfected waveform or a settled limit cycle.

| Maximum step | Closures / reopenings | Accepted substeps | Added numerical loss |
| --- | ---: | ---: | ---: |
| 10 µs | 3 / 3 | 2402 | 21.2596 nJ |
| 5 µs | 3 / 3 | 4800 | 6.5442 nJ |
| 2.5 µs | 3 / 3 | 9600 | 1.6754 nJ |
| 1.25 µs | 3 / 3 | 19200 | 0.4264 nJ |

The three successive scaled endpoint differences are 1.15174, 0.266661, 0.0435791. Their refinement ratios are 4.319, 6.119.
The scaling is stated in the receipt; these are mixed-state norms, not a
percentage error for every component. The remaining near-closure phase and
amplitude differences must be retained as uncertainty. At the finest step,
all six constraint changes complete without subdivision and added numerical
loss is about 0.4264 nJ over the 24 ms window.

See [the trajectory and regression receipt](continuous_phonation_motion.json)
for actual completed intervals, explicit failures, refinement and provenance.
The initial state of each continuation is recorded. Reconstructing that saved
physical state is not a prescribed oscillation: subsequent motion follows the
same gas/tissue equations under constant pressure.

## Validation and remaining work

Ten focused CTests pass in release and address/undefined-sanitizer builds.
They include the two near-steady failures, restart invariance, both free contact/reopening controls,
and existing storage/flow/loss/coupled controls. Fresh default-profile motion
runs compare with historical fixtures, and two 2,000-step contact runs check
that the roundoff fix retains the previous free-motion behavior and work/mass
accounting. The older full contact-refinement study remains historical evidence;
it has not been relabeled as a fresh six-grid run.

Longer nonlinear settling, source/tract coupling, streamed actuator work, body
calibration and new speech are still required. The contact cycles above are
checked transient motion under the stated configuration. The onset results
apply to the symmetric tissue subspace and the explicit full-span shape. They
do not cover every asymmetric tissue mode or validate the local steady BGK
approximation for unsteady separated glottal flow. Performance optimization
remains parked; these are mechanics references.

## Reproduction

Build with CMake, then use a Python environment with NumPy, SciPy and mpmath:

```sh
OPENBLAS_NUM_THREADS=1 python tools/continuous_phonation_reference.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_phonation_modal_experiment.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_phonation_validation.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_phonation_ode.py .06 2e-7 open-ode 32
OPENBLAS_NUM_THREADS=1 python tools/continuous_phonation_ode.py .06 1e-8 open-ode-fine 64
OPENBLAS_NUM_THREADS=1 python tools/continuous_phonation_experiment.py
```

The two physical continuation seeds are recorded in the motion receipt. The
ODE and contact experiment load those seeds if the corresponding local input
files are absent. The motion receipt collector additionally reads the saved
native open-motion, developmental failure, regression and CTest outputs; it is
an audit of those explicit runs, not a replacement for executing them.

# Continuous-source tip contact and reopening

This extends the continuous gas/tissue reference from `dbc880f`. The explicit
full-span linear profile now completes a free closing/reopening motion at each
end, and two controls open from an exactly zero tip at rest. These are source
mechanics controls. They have not produced a new speech recording or changed
the listening assessment of the spatial-pressure alphabet.

## Geometry determines the closing limit

The previous profile has constant-area strips outside the two tissue centres.
A positive lower cap of area `a` therefore remains a parallel strip of length
`c0`, however small `a` becomes. Consider constant enthalpy, relative density
`r`, zero inlet flow and prescribed lower area rate `V`. Strong continuity gives
`Q(z) = -r V z` in that strip. Its exact kinetic energy and BGK loss are

```
Kcap = rho0 r V² c0³ / (6 a)
Pcap = rho0 r vth width V² c0³ / (3 a² G(delta))
delta = (a/width) rho0 r vth / (2 viscosity).
```

`G(delta)` grows logarithmically in the small-gap limit. The positive strip
therefore has divergent kinetic energy and loss under this prescribed motion.
Deleting the strip at zero area would change the admissible motion abruptly.
A finite BGK law at a tapered contact point does not cure this parallel-strip
limit. This argument is conditional on the prescribed constant-enthalpy motion;
it is not a proof that every free compressible trajectory must rebound.

`ContourPressureProfile::linear_span` makes a different, explicit geometry
assumption: `A(z) = max(a(1-z/H) + b z/H, 0)` across the entire thickness `H`.
The two area controls now describe axial ends rather than tissue-centre
plateaus. Tissue masses, springs and contact potential remain the configured
four-mass model. This shape is not measured anatomy or a calibrated vocal fold.
The default remains `flat_caps`.

For lower point contact at `a=0`, with `Qin=0` and the same constant-enthalpy
area rate, `Q(z) = -r V H (t - t²/2)`, where `t=z/H`. Thus

```
K = 11 rho0 r V² H³ / (96 b)
P = rho0 r vth width V² H³ / b²
    * integral_0^1 (1-t/2)² / G(delta(b t)) dt.
```

Both are finite. The native full matrices reproduce the exact kinetic result
and an independently integrated BGK loss result. The fixture uses `H=.003 m`,
`b=3e-6 m²`, `V=-.0003 m²/s`, and enthalpy `700 Pa`. Results and positive/negative
gap approach are in [the profile-limit receipt](contour_profile_limits.json).
The full-span storage, mass derivatives and forces also pass
[11 independent 80-digit cases](linear_contour_storage_reference.json).

## A dissipative update with explicit work accounting

The optional `ContinuousGlottalUpdate::dissipative` update uses the endpoint
shared inverse inertia `G1^-1`, with velocity `v = G1^-1 p1` and displacement
`dq = dt v`. Its configuration discrete gradient accounts for the exact change

```
U1 - U0 + 1/2 p0^T (G1^-1 - G0^-1) p0.
```

The correction is formed in the actual gas-mass and area coordinates, preserving
the cyclic common transferred-volume coordinate. The tissue spring force uses
its exact midpoint discrete gradient. BGK and tissue resistance use the endpoint
configuration. The momentum update then has the work identity

```
H1 - H0 + N + fluid loss + tissue loss + load work
    - reservoir work - external work = residual
N = 1/2 (p1-p0)^T G1^-1 (p1-p0) >= 0.
```

`N` is numerical damping and is reported separately from physical dissipation.
Conservation is insufficient evidence of accuracy; the refinement controls
check motion and decreasing `N`. The symmetric update remains available for
fixed topology. A flat-cap topology change is rejected even in dissipative
mode. Neither mode removes gas mass or inserts a source waveform.

At a closed exterior port the corresponding transferred-volume displacement
and velocity are eliminated exactly before evaluating storage. An earlier
version only enforced a Newton residual tolerance, allowing tiny nonzero
average closed-port flows; the strict control caught that error. At exact
contact and rest, the initial iterate uses a pressure/spring acceleration
predictor, not a positive-gap floor. The exact-contact kinetic derivative has
a defined closed-side value. Small new ports use a Jacobian perturbation bounded
by the acoustic volume their actual area can transport during the time step.

## Resolving the physical interval without a gap floor

Positive taper intervals are partitioned by factors of four in area. This
resolves the full thin-gap region; no interval is discarded. For the full-span
profile, integration runs from the narrow end. The coordinate reflection
`z -> H-z` swaps pressure moments and wall rates, and reverses the three flux
coordinates. The matrices undergo that same signed permutation.

This addresses a failure found when opening the upper tip from exact rest:
subtracting nearly equal axial coordinates destroyed the narrow-end quadrature
accuracy. Lower-tip integration already used a coordinate beginning at zero.
Both orientations now use that representation, including loss integration.
Reversal checks cover positive areas down to `1e-20 m²`, zero, and negative
contact coordinates. The tiny-area checks establish numerical representation,
not physical validity of this continuum/BGK approximation at those dimensions.
Exactly zero endpoints are also assigned their exact basis coordinates rather
than reconstructed by a cancellation-prone root expression.

## Validation and remaining work

The native free-motion controls start one area at `1e-7 m²`, the other at
`3e-6 m²`, and the closing tissue pair at `-.05 m/s`. The configured rest area
is `3e-6 m²`; reservoir and downstream efforts are zero. Both controls retain
positive gas moments, seal one end, and reopen it over two milliseconds.
Refinement uses 2,000, 4,000 and 8,000 steps. The receipt below records actual
closure/reopening times, final-state differences, physical and numerical loss,
and mass/work residuals. It does not impose a claimed asymptotic order at the
singular contact transition.

| Closing end | Finest closure | Finest reopening | Final-state refinement ratio | Added loss, 2k / 4k / 8k steps |
| --- | ---: | ---: | ---: | ---: |
| lower | 78.50 µs | 809.50 µs | 1.859 | 2.518 / 1.460 / 0.836 nJ |
| upper | 91.00 µs | 436.25 µs | 1.928 | 1.400 / 0.890 / 0.545 nJ |

Across all six runs the maximum cumulative work residual is 6.15e-17 J;
the maximum sampled mass-balance residual is 1.54e-21 m³ in reference-volume
mass units. The last column isolates solver damping; it is not BGK or tissue loss.

Separate exact-zero tests start with zero tissue velocity and `700 Pa` internal
enthalpy, a `700 Pa` upstream reservoir and zero downstream effort. Each opens
its initially closed tip in a `1e-7 s` step without a prescribed opening.
Seven focused release tests and seven address/undefined-sanitizer tests cover
storage, flow, loss, coupled equilibrium and the two contact controls. Four
fresh default-profile motion runs and a constitutive run are compared against
saved earlier fixtures; those regressions reuse the earlier oracle evidence.

See [the contact/refinement and validation receipt](continuous_contact_reference.json).
Reproduction tools are `continuous_contact_experiment.py`,
`continuous_contact_reference.py`, `continuous_contact_regression.py`,
`contour_profile_limit_reference.py`, and `linear_contour_storage_reference.py`.
The native controls are registered with CTest. To reproduce the new references
from a configured release build:

```sh
cmake --build build -j4
OPENBLAS_NUM_THREADS=1 python tools/contour_profile_limit_reference.py
OPENBLAS_NUM_THREADS=1 python tools/linear_contour_storage_reference.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_contact_experiment.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_contact_regression.py
```

Python requires NumPy, SciPy and mpmath. The receipt collector additionally
requires the saved release/sanitizer CTest logs; it checks their seven-test
completion and hashes the exact source, binaries and output files.

Sustained pressure-driven oscillation, a state with the whole passage closed,
moving rest-area actuator work, source/tract coupling, and new speech remain
unfinished. This expensive quadrature/Newton implementation is a mechanics
reference, not a realtime speech implementation. Performance optimization
remains parked by user instruction. The local steady diffuse-wall BGK closure,
full-span profile and four-mass tissue approximation still need physical
validation; successful numerical controls do not establish human realism.

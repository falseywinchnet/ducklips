# Radiation at the continuous-source airway outlet

The developing continuous gas/tissue source and moving mouth now have an
explicit passive radiation boundary. It loads the outlet flow, stores near-field
air energy and accounts separately for energy radiated away and work from a
changing aperture. This is an opt-in research connection, not a new speech
recording or a completed utterance renderer.

## Boundary model and energy

`RadiationBoundary` uses the existing area-equivalent monopole/end-correction
approximation. For aperture area `A`, air density `rho`, sound speed `c`, end
correction `e = 0.7` and radiation factor `b = 1.4`,

```
M = rho e sqrt(A/pi) / A
H = 1/M
K = b rho / (4 pi c)
R = M^2/K
Z(s) = M s / (1 + s M/R)
```

The low-frequency reactance is `M omega`; the resistance approaches `K omega^2`.
These are declared model parameters. This is not a calibrated high-frequency
solution for the actual lip shape.

For duct canonical momentum `p` and reciprocal inertia `G`, duct flow is `G p`.
Radiation momentum `r` has near-field energy `H r^2/2`. The pressure drop and
fixed-geometry exchange are

```
P = R (G p - H r)
p' = -P
r' = P
radiated power = P^2/R >= 0
```

The exchange conserves `p+r`. Its difference mode decays at `R(G+H)`, so it is
advanced exactly with an exponential, including its integrated radiated energy.
`expm1` retains small increments without subtracting nearly equal numbers.
There is no explicit small-timestep stability restriction for this isolated
exchange; coupled trajectory accuracy still requires refinement.

Positive-area changes retain canonical radiation momentum. The resulting
`(H_new-H_old) r^2/2` is geometry work, not radiation loss. At exactly zero area,
`H` and near-field energy are zero and the absent radiation memory is retired.
The closed boundary accepts only zero duct conductance and radiation momentum;
reopening starts without a ghost exterior state. The adapter continues to own
the separate retirement of closed duct momentum. No minimum aperture is added.

## Connection and continuation

`GraphContinuousGlottis::advance` accepts an optional `RadiationGeometry` with
the exterior branch, its endpoint conductances, and actual midpoint and endpoint
aperture areas. The boundary uses a symmetric split:

1. Change its metric to the midpoint aperture and exchange momentum for half a step.
2. Advance the existing implicitly coupled source and graph for the full step.
3. Exchange for another half step and change to the endpoint aperture.

The first and last exchanges use their respective endpoint duct metrics.
Radiation is not part of the implicit middle solve. The independent continuous
ODE below measures this splitting error. Prepared advances allocate no memory.
A failed trial preserves the accepted source, graph, radiation and receipts.

The planar adapter enables the connection with

```sh
ducklips_planar_inertive_motion_probe --lip-profile --epilarynx \
    --continuous-folds bodies/reference.body --radiation
```

The aperture is the integrated lip profile's exterior area, not its last cell's
midpoint area. The existing duct half-section remains physical interior air;
the radiation end correction is represented only by the new boundary state.
Continuation appends radiation momentum after the source's thirteen physical
values. Its area, inverse inertia, momentum and energy round-trip exactly in
the tested interval continuations. Source reconstruction retains its separately
measured floating-point round-trip error.

The middle graph's average flow is the step's volume transfer, including at the
source port. Radiation exchanges alter momentum without separately transferring
nodal air. Combined stored-air accounting therefore retains the same reservoir
inflow and exterior lip outflow as before.

## Independent controls

The [kernel reference](radiation_boundary_reference.json) compares 100 seeded
cases against a separate matrix exponential in energy-scaled coordinates.
Maximum normalized discrepancies are 6.76e-15 in state, 1.45e-14 in radiated
energy and 3.86e-16 in prescribed aperture work. Native controls cover another
48 combinations, exact closure, reopening, invalid inputs and failure rollback.

The [coupled reference](continuous_radiation_reference.json) solves independent
time-continuous equations for the open gas/tissue source, a fixed cyclic
three-node tract and either constant or smoothly changing aperture. A Radau
solve with 48-point source quadrature and tighter tolerance differs from the
32-point reference by less than 1.82e-9 in the declared normalization. The native
two-millisecond trajectories use 200, 400, 800 and 1,600 steps.

| Finest run / last halving | Constant aperture | Changing aperture |
|---|---:|---:|
| Endpoint scaled error | 1.343e-4 | 1.179e-4 |
| Endpoint error reduction factor | 3.991 | 3.992 |
| Radiated-energy error | 1.512e-9 J | 1.344e-9 J |
| Radiated-energy error reduction factor | 3.990 | 3.991 |
| Geometry-work error | 0 J | 7.227e-12 J |
| Maximum sampled scaled state error | 0.03712 | 0.03712 |

The startup transient remains underresolved even at the finest 1.25 microsecond
step. Small endpoint errors do not establish accuracy throughout the trajectory.
For this fixed graph, source plus graph plus near-field energy closes against
source losses, radiated energy, external work and aperture work: the maximum
residual over all runs is 6.34e-15 J. This is not a claim of a complete energy
identity for the separately moving planar graph.

## Actual moving-mouth controls

The [moving radiation receipt](planar_radiation_reference.json) covers the same
separately initialized Bravo closing and opening controls as the preceding
[moving-airway work](CONTINUOUS_MOVING_AIRWAY.md). Each continues state across
two one-millisecond intervals. Twelve runs at 256, 512 and 1,024 subdivisions
per millisecond total 7,171 accepted steps.

Independent Boolean geometry, lip quadrature, graph stage solves and radiation
exchanges check 3,079 geometry snapshots at the first refinement. The graph
reference uses the recorded coupled source flux. This validates geometry and
coupling; it is not an independent whole-source/moving-mouth ODE. Finer runs
check ledgers and temporal differences.

| Measured maximum | Value |
|---|---:|
| Independent graph pressure discrepancy | 2.00e-9 Pa |
| Independent radiation momentum discrepancy | 3.55e-14 Pa s |
| Independent step radiated-energy discrepancy | 1.84e-19 J |
| Independent step radiation-geometry-work discrepancy | 6.62e-20 J |
| Local combined stored-air residual | 9.45e-21 m3 |
| Per-interval combined stored-air residual | 7.59e-20 m3 |
| Source step work residual | 4.54e-17 J |
| Source/tract port-work mismatch | 5.30e-23 J |

Fully sealed lip intervals have exactly zero exterior flow and radiation loss.
A step ending at closure may still transport air before that endpoint.

| Refinement, subdivisions/ms | Closing endpoint change | Opening endpoint change |
|---|---:|---:|
| 256 to 512 | 51.96 Pa | 8.09 Pa |
| 512 to 1,024 | 39.54 Pa | 2.76 Pa |

Neither control is established as temporally converged. Radiation barely
changes closing sensitivity relative to the preceding ideal-outlet experiment.
At closure the cavity cannot radiate through the lips; interior viscous,
thermal and yielding-wall losses remain absent in this path. Their addition
must preserve the physical meaning of stored air, including any thermal or
wall-displacement states. Generic pressure damping would not meet that condition.
General changing-region transport, physiological actuation and whole-utterance
integration also remain unfinished. No intelligibility, realism or realtime
performance claim follows from these controls.

## Validation and reproduction

Thirteen focused CTests pass in Release and thirteen under ASAN/UBSAN. All four
actual moving radiation intervals pass under ASAN/UBSAN at 256 subdivisions/ms.
Eleven historical adapter outputs reproduce their entire JSON exactly: seven
older geometry/airway cases and four non-radiating continuous-source intervals.
The [validation receipt](continuous_radiation_validation.json) records fresh
test logs and verifies the nested input, output, binary and source hashes.

```sh
cmake --build build --target ducklips_radiation_boundary_probe \
    ducklips_graph_continuous_glottis_probe ducklips_planar_inertive_motion_probe
python tools/radiation_boundary_reference.py
python tools/continuous_radiation_reference.py
python tools/planar_continuous_reference.py --radiation
python tools/radiation_regressions.py
```

Use the scientific Python environment with NumPy, SciPy and Shapely. The planar
checks use the saved native pose/profile inputs pinned in the receipt. Radiation
outputs have separate directories and do not overwrite the previous non-radiating
controls or the spatial-pressure listening recording.

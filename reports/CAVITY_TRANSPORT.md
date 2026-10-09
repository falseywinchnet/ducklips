# Conserving air when a cell appears or disappears

The native mean-pressure solver now permits a disappearing cell to expel its
stored air without retaining artificial volume or compliance. It also fills
a newly exposed cell through flow instead of copying a neighbour's pressure.
Analytical and independent numerical checks pass for the conservation equation.

The complete Foxtrot transition is **not ready for speech integration**.
A 66-pose audit exposes split cells missed by the earlier four snapshots, and
independent interpolation of endpoint volume and aperture has not established
pressure convergence near tiny cells. These are retained failures, not passed
realism gates. No new speech recording or enhanced spectrogram is attributed
to this work; the latest [listening feedback](LISTENING_SPATIAL_2026_10_09.md)
still reports November improved and Foxtrot jammed around /tr/.

A subsequent [native region decomposition](PLANAR_REGIONS.md) resolves every
split-cell pose in this audit and preserves the /t/ seal in a fixed-pressure
control. The temporal refinement failure below remains unresolved; the
combined-cell experiments are retained as failure evidence.

## Zero storage without deleting air

Each persistent Eulerian cell uses the existing linearized storage variable
`M = V + C p`, in reference-volume units. The step equation is

```
M_new - M_old + dt * sum(outgoing interval flows) = dt * inflow.
```

For a completely occupied cell, `V_new = C_new = 0`, so its old contents must
leave through the passages available during the step. Its pressure then acts
as an algebraic transport multiplier; it is not a physical endpoint pressure
or stored acoustic energy. No minimum volume or epsilon compliance is added.
An already absent, isolated cell with no inherited mass or source has its
unused pressure fixed to zero. A cell cannot disappear with trapped nonzero
mass. The failed step leaves the accepted pressure, volume, compliance and
flow state unchanged.

Every component of interval-open passages must have positive storage or an
open pressure reservoir. The pressure solver checks this before Newton.
The spanning-forest preconditioner also handles a zero-storage cell whose
only remaining exit is a non-tree edge, without changing the physical graph.
Existing positive-compliance speech equations retain their previous route.

This API does not infer interval openings from endpoint geometry. Giving a
vanishing cell only its already-closed endpoint faces cannot expel its air.
The caller must supply a consistent geometry/flow law over the step. The
current research drivers test that distinction explicitly.

## Analytical checks

A two-cell disappearance test fixes the expelled flow independently from
`Q = M_old / dt`. The remaining chamber receives that same stored amount.
The nonlinear variant also matches the independently calculated viscous and
quadratic pressure drop, with pressure error 1.16e−8 Pa. Reappearing-cell
pressure agrees with a separate two-equation calculation. Sealed disappearance
is rejected and the original state remains usable. Stepping allocates nothing.

A moving wedge provides a geometry-consistent birth/death case. Let its
volume be `V = a h²` and conductance `G = b h`, where `h` is elapsed or
remaining time (wedge height is proportional to it), and `C = V/B`, with
`B = rho c²`. With a zero-pressure reservoir,
constant pressure solves the equations exactly:

```
opening: p = -2a / (b + 2a/B)
closing: p =  2a / (b - 2a/B).
```

Using endpoint volumes and the interval-average conductance preserves this
solution. At 1, 2, 4, 8 and 16 steps, the native fixture has zero observed
pressure error and maximum reference-volume balance error 1.66e−24 m³.
This establishes a useful finite limiting case at cell creation/destruction;
it does not validate the tongue's interpolated geometry.

## Full geometry audit: a real false connection

The native geometry probe can now export each selected pose rather than only
the last one. Shapely/GEOS independently checks all 66 poses from 5.020 through
5.085 s, including polygon validity, every cell volume, split fragments and
whole-domain connectivity. All supplied polygons are valid; no automatic
geometry repair is used. Across this interval there are 195 cell-appearance events and 176
disappearance events.

Eight poses fail the one-pressure-per-cell requirement:

| Time | Cells containing disconnected fluid pieces |
|---|---|
| 5.032 s | 3760 |
| 5.037 s | 3630 |
| 5.079 s | 2853, 2917, 2981, 3045 |
| 5.080 s | 2853, 2917, 2981 |
| 5.081–5.083 s | 2853, 2917 |
| 5.084 s | 2853 |

At 5.032 s, cell 3760 contains pieces belonging to different sealed chambers.
Their volumes are 0.000351780 and 0.000004694 mL. The exact fluid domain has
two components, while the graph that combines these pieces into one cell has
only one. This would short-circuit the /t/ closure. The other split-cell
cases can also create incorrect local routes even when the whole-domain
component count happens to match. The driver rejects every split pose.

## Partial transport replay and failed refinement

The two audited intervals 5.020–5.031 and 5.038–5.078 s are tested separately,
starting at uniform zero pressure, with a sealed outer boundary. They retain
all 4,096 persistent grid locations and all 8,064 possible internal edges.
No pressure is carried from an excluded pose into these tests. They are not
a replay of the pressure state in the assessed speech recording.

The numerical fixture specifies `G = mobility * aperture / distance`, with
mobility 1e−7 m²/(Pa s). It does not pretend that computational faces are tube
walls or that this mobility is calibrated vocal-tract friction. Endpoint
volumes and apertures are interpolated independently, and each step uses a
trapezoid-average aperture. Runs use 1, 2 and 4 substeps per saved millisecond.

| Interval | Appearing / disappearing cells | Max native/reference pressure difference | Max total reference-volume discrepancy |
|---|---:|---:|---:|
| Before closure | 17 / 15 | 2.26e−7 Pa | 5.43e−19 m³ |
| Closure through release | 147 / 136 | 5.39e−6 Pa | 3.62e−17 m³ |

The independent sparse solve evolves its own state rather than restarting
from the native solution. These results establish numerical agreement and
storage conservation for the supplied equations, not valid speech physics.

Refinement gives mixed evidence. Volume-weighted RMS endpoint differences
fall from 81.8 to 42.2 Pa in the first interval and 591.8 to 265.9 Pa in the
second. However, maximum differences grow from 2,515 to 3,534 Pa and from
6,585 to 7,561 Pa respectively. The worst cells are tiny: 3.01e−12 and
2.70e−11 m³. Model density ratios across the experiments range from 0.279 to
2.827, far outside the small perturbations assumed by this linearized model.
Pressure refinement and physical adequacy are therefore unproven.

A geometric reason to investigate is that independently linear volume and
aperture interpolation need not preserve the relationship between swept
volume and an opening. For example, near the birth of a wedge, volume is
quadratic in height and aperture is linear. Replacing both with unrelated
linear interpolants changes the limiting filling problem. The exact wedge
test verifies the solver when those quantities are consistent; it does not
prove that this explains every Foxtrot pressure excursion. Reconstructing
substep boundaries and their swept openings is required to resolve that.

## Evidence and remaining work

The full Release suite passed 14/14, and the affected stream, cavity-network
and mean-airway tests passed ASan/UBSan, 3/3. The subsequently expanded wedge
and nonlinear-disappearance fixtures also passed directly in Release and
under ASan/UBSan. Both unrefined trajectory inputs passed sanitizers; their
maximum pressure difference from Release was 3.23e−6 Pa. These are numerical
and memory-safety checks, not acceptance of the failed geometry/refinement.

[cavity_transport_reference.json](cavity_transport_reference.json) includes
the geometry audit, rejected poses, refinement results, commands and hashes.
The tool intentionally produces a diagnostic report with
`whole_transition_valid: false` and `pressure_refinement_established: false`.
The existing fixed-graph audit was rerun against the extended implementation.
Local detailed outputs are under `out/cavity-transport/`.

Component-resolved spatial geometry is now independently checked at all 66
poses. The next work is consistent geometry and conservative region-state
correspondence over the time interval. The mean-pressure model then still needs physical
flow laws, acoustic state and moving-boundary work, preserved nasal/side
cavities and glottal/radiation coupling, and new listening evidence. The blade
and transverse anatomy also remain unvalidated. No speech-realism or speed
gate has been relaxed or declared satisfied.

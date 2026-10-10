# Heat diffusion with conserved air in the continuous-source tract

The developing source/tract path now includes explicit finite-gap thermal
storage. It changes pressure as air exchanges heat with isothermal walls while
retaining the additional mass-storage term. It is connected to the continuous
gas/tissue source, radiation boundary, and eight throat sections in a declared
parallel-wall research configuration. This is not a new speech rendering.

## Why pressure damping alone is wrong

Let `C = V/(rho c^2)` be adiabatic acoustic compliance, `P` pressure, and `gamma`
the heat-capacity ratio. The graph's primary charge is `q = C P`, in reference
cubic metres of air. After thermal equilibration, the same pressure stores
`gamma C P`. The difference must remain part of the conserved air, rather than
being erased from the state when pressure relaxes.

The new representation assigns thermal charges `z_j` to passive diffusion
modes. With positive capacity ratios `kappa_j` and decay rates `lambda_j`,

```
z_j' = lambda_j (kappa_j q - z_j)
q'   = -sum(z_j')                 [during local thermal exchange]
total stored acoustic air = q + sum(z_j)
energy = q^2/(2 C) + sum[z_j^2/(2 C kappa_j)]
thermal dissipation = sum[lambda_j (kappa_j q-z_j)^2/(C kappa_j)] >= 0
```

This is acoustic energy dissipation; it is not an absolute calorimetric heat
measurement of the body. The thermal charge can have either sign. No nodal air
flows through a wall in this exchange.

For a slit of half-gap `a`, the exact constant-wall-temperature diffusion
response is

```
Upsilon(s) = tanh(sqrt(s a^2/D)) / sqrt(s a^2/D)
C_effective(s) = C [1 + (gamma-1) Upsilon(s)]
```

Here `D` is the existing thermal diffusivity. This response follows the linear
heat equation across the slit. Its parallel-wall, plane-wave and slowly varying
cross-section assumptions are the narrow-region assumptions described in the
[COMSOL model documentation](https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.183.html).
The present implementation does not solve the full moving-wall heat/advection
PDE or calibrate a human throat's cross section.

The exact thermal spectrum has weights `w_n = 2/[pi^2(n+1/2)^2]` and
dimensionless poles `p_n = pi^2(n+1/2)^2`, with `sum(w_n)=1`.
`kappa_j=(gamma-1)w_j` and `lambda_j=p_j D/a^2`. The first eight modes are exact;
the rest are grouped into positive intervals, including the infinite tail.
Each group preserves its total weight and first low-frequency moment:

```
w_group = sum(w_n)
p_group = sum(w_n) / sum(w_n/p_n)
```

The 24 groups therefore retain the exact static compliance, rather than a
low-frequency cutoff chosen independently of the gap. The coefficient table is
independently regenerated using polygamma sums in the reference tool.

## Native step and geometry

`SlitThermal` analytically eliminates the conserved total charge. At fixed
geometry the remaining symmetric positive-definite generator is

```
S = diag(p) + sqrt(kappa*p) sqrt(kappa*p)^T
```

Its modes are prepared once. An advance evaluates their exponential decay and
integrated dissipation. `expm1` preserves small changes. The zero decay mode is
handled through the conserved total rather than a nearly zero numerical
eigenvalue. Prepared advances allocate no memory.

Thermal charges persist through a geometry change. Their compliance metric
changes with the node volume, and the resulting stored-energy change is
`thermal_geometry_work`. A changing gap alters the diffusion rates. Empty or
nonpositive storage is rejected: disappearing regions require a separate
physical transport operation. There is no minimum-volume or minimum-gap floor.

`GraphContinuousGlottis::prepare_thermal` declares the affected nodes and their
initial compliance/gap. Explicit equilibration initializes an isothermal
pressure state; zero thermal charge represents an initially adiabatic pressure
disturbance. The graph wrapper applies a thermal half step before the existing
source/graph advance, then changes the thermal endpoint metric and applies its
second half step. Thermal and radiation exchanges act on different primary
coordinates in those exterior half steps. They do not replace the source's
implicit pressure/flow coupling.

The graph's reported average flow remains its middle-step mass transfer.
Conservation must include both graph charge and thermal charge. Trial failure
preserves source, graph, radiation, all thermal memories and receipts. Late
invalid thermal geometry is checked after the source trial to exercise rollback.

## What the independent checks show

The [local diffusion reference](slit_thermal_reference.json) compares the
grouped transfer function with the exact `tanh` expression over dimensionless
angular frequencies from 1e-6 through 1e8. Maximum relative error in the full
complex acoustic compliance is **7.975e-5**. This is the error in the full
compliance, not relative error in its smaller thermal contribution.

One hundred random states are also advanced with a separate full 25-coordinate
matrix exponential. Maximum normalized state discrepancy is 5.40e-14 and
relative energy discrepancy 3.81e-14. Twenty-seven native controls exercise
widely separated compliances, gaps and times. Their largest relative energy
residual is 1.44e-15. A closed-cavity control starting at 700 Pa with zero thermal
charge reaches 500 Pa for `gamma=1.4`, retaining total air. Additional controls
check geometry work, equilibrium and invalid-input rollback.

The source/graph/radiation/thermal control is checked against an
[independent continuous ODE](continuous_thermal_reference.json), including the
thermal state and integrated dissipation. Its 48-point/tighter-tolerance source
reference differs from a 32-point control by 3.72e-9 in the declared state
normalization. The run is 0.5 ms, open-source, with a fixed cyclic three-node
tract, fixed thermal gaps and fixed aperture radiation.

| Steps | Endpoint scaled state error | Maximum sampled scaled state error | Thermal dissipation error |
|---|---:|---:|---:|
| 200 | 2.237e-3 | 5.969e-2 | 2.563e-9 J |
| 400 | 5.611e-4 | 3.710e-2 | 9.421e-10 J |
| 800 | 1.954e-4 | 1.481e-2 | 3.420e-10 J |
| 1,600 | 1.857e-4 | 3.996e-3 | 1.232e-10 J |
| 3,200 | 8.902e-5 | 1.004e-3 | 4.409e-11 J |

The [additional 3,200-step check](continuous_thermal_refinement.json) is a fresh
native run evaluated at the exact sample times of the saved independent ODE.
It reuses that verified reference; it does not claim a newly integrated ODE.
Independent thermal dissipation is 7.262476446e-7 J. All five runs retain
combined work residual below 6.06e-16 J and stored-air residual below 4.62e-21 m3.

The maximum sampled error falls under refinement. The endpoint error does not
uniformly decrease by four, and these runs do not establish an asymptotic order
for the split thermal/source trajectory. At 1,600 steps the largest endpoint
normalized error is in source pressure, rather than a thermal memory. This
accuracy limit remains visible; an energy ledger cannot substitute for it.

## Actual moving-airway connection

The planar adapter explicitly enables thermal storage with

```sh
ducklips_planar_inertive_motion_probe --lip-profile --epilarynx \
    --continuous-folds bodies/reference.body --radiation --thermal-slit-width 0.03
```

This configuration treats the eight existing throat sections as slits of the
declared 30 mm width, using `half_gap = area/(2 width)` and their existing volumes.
It preserves the axial geometry and source attachment. The width is a research
configuration, not a measured anatomical fact. Only those throat nodes receive
thermal storage; mouth/lip thermal geometry, side-wall heat transfer, viscous
drag and yielding walls remain absent in this path.

Initial thermal charge equilibrates with each declared initial throat pressure.
Continuation appends 24 thermal charges per site after the source and optional
radiation state. All thermal state and geometry fields round-trip exactly at
the tested interval boundaries; source reconstruction retains its
separately measured floating-point round-trip behavior.

The [moving control receipt](planar_thermal_reference.json) covers twelve runs
of the same separately initialized Bravo closing/reopening intervals. It checks
7,171 accepted steps and 3,079 independent geometry snapshots. The reference
uses independent Boolean geometry, lip quadrature, a graph stage solve and
full thermal matrix exponentials, driven by the recorded source flux. It is
not an independent ODE of the entire moving source/mouth system.

| Measured maximum | Value |
|---|---:|
| Independent graph pressure discrepancy | 7.80e-9 Pa |
| Independent thermal charge discrepancy | 5.13e-20 m3 |
| Independent step thermal-dissipation discrepancy | 3.39e-18 J |
| Local combined stored-air residual | 8.98e-21 m3 |
| Per-interval combined stored-air residual | 5.70e-20 m3 |
| Source step work residual | 5.22e-17 J |
| Source/tract port-work mismatch | 5.30e-23 J |

The Python audit originally inferred an integer dtype from initially zero
thermal charges, discarding small nonzero values after the first step. The
reference now explicitly uses floating-point storage; the native thermal kernel
did not require a change for that audit failure.

Closing endpoint differences remain 51.96 Pa then 39.54 Pa under step halving;
opening differences are 8.09 Pa then 2.76 Pa. The throat thermal path has little
effect on those refinement differences. At the finest resolution its accumulated
dissipation is 4.53e-8 J over the two closing intervals and 2.45e-8 J over the two
opening intervals. Throat volumes are fixed in these particular controls, so
their thermal geometry work is zero; changing-compliance work is covered by the
separate native kernel control.

## Validation and remaining work

Fifteen focused CTests pass in Release and fifteen under ASAN/UBSAN. All four
actual moving thermal/radiation intervals pass under ASAN/UBSAN at 256
subdivisions/ms. Fifteen historical JSON outputs reproduce exactly, including
the previous radiation path. The [validation receipt](continuous_thermal_validation.json)
verifies nested source, input, output, binary and test-log hashes.

```sh
cmake --build build --target ducklips_slit_thermal_probe \
    ducklips_graph_continuous_glottis_probe ducklips_planar_inertive_motion_probe
python tools/slit_thermal_reference.py
python tools/continuous_thermal_reference.py
python tools/planar_continuous_reference.py --radiation --thermal-slit-width 0.03
python tools/radiation_regressions.py --output out/thermal-regressions --include-radiation
```

Use the existing scientific Python environment with NumPy, SciPy and Shapely.
The stored native input in the additional refinement receipt can be passed to
`ducklips_graph_continuous_glottis_probe --thermal-fixture` to reproduce it.

This implements a physical heat-diffusion storage path, not all interior losses.
Viscous wall drag is the next missing dissipative connection. The split coupling's
time accuracy, broader thermal geometry, region transport, physiological
actuation and whole-utterance integration still need work. There is no new
intelligibility, realism or realtime-performance claim, and the spatial-pressure
listening recording remains unchanged.

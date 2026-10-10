# Finite-gap viscous drag sharing the tract's existing flow

The developing continuous-source path now includes stationary-wall slit Stokes
drag and transverse velocity-profile memory. Two physical half-sections can
load the same graph flow without duplicating its plug-flow inertia. The model
is connected to the thermal/radiation/source graph and the seven internal throat
branches of the actual moving-mouth controls. It is not a new speech recording.

## Physical impedance and state

For a uniform slit segment with plug-flow inertia `L`, half-gap `a`, kinematic
viscosity `nu`, and Laplace variable `s`, integrating the linear no-slip Stokes
velocity profile across the slit gives

```
Z(s) = L s / [1 - tanh(sqrt(s a^2/nu))/sqrt(s a^2/nu)]
```

This is a stationary parallel-wall, linear plane-wave approximation. Its
cross-section and wavelength assumptions have the same narrow-region scope as
the [thermoviscous duct models documented by COMSOL](https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.183.html).
It does not include axial wall motion, finite-amplitude convective flow or
anatomical calibration.

Writing `t = s a^2/nu`, the dimensionless impedance is

```
t/[1-tanh(sqrt(t))/sqrt(t)] = t + 3 + sum[2t/(t+x_n^2)]
tan(x_n) = x_n,    n pi < x_n < (n+1/2) pi,    n=1,2,...
```

The graph already owns `L s`. The added model supplies only

```
R0 = 3 L nu/a^2
R_j s/(s+lambda_j), with R_j = L nu/a^2 w_j, lambda_j = nu/a^2 p_j
```

The steady resistance is the slit Poiseuille resistance. The low-frequency
profile contributes a further `L/5` of effective inertia, so a steady parabolic
profile has 1.2 times the plug-flow kinetic energy at the same volume flow.
That extra energy belongs to profile memory; the ordinary air-column inertia
is not counted twice.

For full-branch reciprocal inertia `G`, graph momentum `p`, flow `Q=G p`, and
segment memory momenta `z_j`, the viscous part of the motion is

```
z_j' = R_j Q - lambda_j z_j
p'   = -R0 Q - sum(z_j')
profile energy = sum[lambda_j z_j^2/(2 R_j)]
dissipated power = R0 Q^2 + sum[(R_j Q-lambda_j z_j)^2/R_j] >= 0
```

These states carry velocity-profile energy, not additional stored air. Combined
air conservation still uses the source gas, actual tract volume, acoustic charge
and the previously added thermal charge.

## Reduced spectrum and exact local decay

`SlitViscous` uses 24 positive modes. The first eight roots remain separate;
subsequent roots are grouped geometrically, with an infinite final tail. For
each group, define `m1=sum(2/x_n^2)` and `m2=sum(2/x_n^4)`; its positive pole and
weight are `p=m1/m2` and `w=m1^2/m2`. This preserves two low-frequency moments.
The exact sums `sum(2/x_n^2)=1/5` and `sum(2/x_n^4)=1/175` retain the entire
tail. The reference tool regenerates the table with 40-digit root calculations.

This reduction has a measured approximation error. Over dimensionless angular
frequencies 1e-6 through 1e8, comparison with the exact slit expression gives:

| Maximum relative error | Value |
|---|---:|
| Full complex impedance | 1.5066e-4 |
| Resistance alone | 1.9387% |
| Viscous correction after subtracting plug inertia | 1.9586% |

The small full-impedance error must not be confused with the larger error in
its dissipative part. This is a bounded numerical approximation to the declared
slit model, not a claim of a perfected loss model or a human throat measurement.

The prepared local system is symmetric positive definite in energy-scaled
coordinates. Its 25 modes include the primary branch momentum. A step evaluates
their exponential decay and integrated dissipation. Small changes use `expm1`;
when all modes decay appreciably, direct modal reconstruction avoids leaving
subtraction residue after long relaxation. There is no timestep clamp or
added global damping.

Changing segment inertia retains its memory momenta and records the resulting
profile-energy change as geometry work. The graph owns work from its own
changing inertia. A changed segment/branch inertia ratio rebuilds the local
decay basis without allocation. Nonpositive paths are rejected; disappearance
of an air path requires explicit physical state retirement or transport and is
not implemented by inserting a minimum area.

## Shared-branch source coupling

`GraphContinuousGlottis::prepare_viscous` declares segment owners and initial
metrics. Segments on a branch must agree on its conductance, and their inertances
may sum to no more than the graph branch's existing inertia. This prevents
assigning overlapping physical lengths silently.

The wrapper applies viscous half steps in forward segment order, advances the
existing source/graph/radiation connection, then applies the ending half steps
in reverse order. Shared-flow segment updates do not commute, so reversing the
order is necessary for the symmetric composition. Thermal half steps surround
this operation as before. The source still solves its implicit pressure/flow
connection; viscous and radiation exchanges are outside that middle solve.

Prepared advances allocate no memory. Failure preserves all accepted source,
graph, radiation, thermal, and viscous states and receipts, including rejection
of late invalid viscous geometry after trial source evolution. Controls also
reject an attempted segment partition exceeding the branch's physical inertia.

## Independent evidence

The [local reference](slit_viscous_reference.json) checks the exact slit transfer
function and 100 random states against independent full matrix exponentials.
Maximum scaled state discrepancy is 1.91e-13 and relative energy discrepancy
2.69e-14. Eighty-one native controls cover different segment fractions, inertias,
gaps and timesteps. Their largest relative energy residual is 3.19e-15. Separate
controls verify profile inertia, geometry work, rollback and complete decay.

The [fully coupled continuous reference](continuous_viscous_reference.json)
integrates an open source, fixed cyclic three-node tract, aperture radiation,
thermal storage and three viscous segments. Two segments share one branch; the
third shares the exterior branch with radiation. The independent ODE evolves
all memory states and integrates viscous and thermal dissipation explicitly.
Its refined quadrature/tolerance references differ by 5.49e-9 in the declared
state normalization.

| Steps over 0.5 ms | Endpoint scaled state error | Maximum sampled scaled state error | Viscous dissipation error |
|---|---:|---:|---:|
| 200 | 2.388e-3 | 5.967e-2 | 2.117e-9 J |
| 400 | 6.007e-4 | 3.709e-2 | 7.930e-10 J |
| 800 | 1.981e-4 | 1.480e-2 | 2.915e-10 J |
| 1,600 | 1.755e-4 | 3.995e-3 | 1.059e-10 J |

Independent accumulated viscous dissipation is 2.903421745e-6 J. The maximum
combined work residual over these runs is 4.19e-16 J and stored-air residual
4.50e-21 m3. Errors decrease, but the last endpoint reduction factor is only
1.129. The split coupled system is not established as asymptotically converged;
passing the energy ledger is not evidence of exact pressure motion.

## Actual throat geometry and scope

The explicit planar research configuration is

```sh
ducklips_planar_inertive_motion_probe --lip-profile --epilarynx \
    --continuous-folds bodies/reference.body --radiation \
    --thermal-slit-width 0.03 --viscous-slit-width 0.03
```

Each of the seven internal throat branches consists of its two actual adjoining
half-sections. The adapter assigns each segment `L=rho length/(2 area)` and
`half_gap=area/(2 declared_width)`, using the branch's existing full conductance.
Thus fourteen profile memories share seven flows. The 30 mm width is an explicit
research configuration; side walls and a measured three-dimensional throat
shape are not represented by it.

Initial profile memories are zero, matching the initialized zero flow.
Continuation appends 24 memory momenta per segment after source, optional
radiation and thermal states. All viscous state and geometry fields retain exact
continuation at the tested interval boundaries.

This addition covers only the internal throat branches. Mouth/lip drag, the
throat's mouth-attachment half-section, the source-owned glottis-to-throat
half-section, yielding walls and disappearing viscous paths remain unfinished.
In particular, the source-owned half-section shares its inertia with gas and
tissue coordinates; it cannot be supplied by an unrelated graph flow state.

The [moving-airway receipt](planar_viscous_reference.json) records the actual
Bravo closing/reopening controls. Their independent graph audit uses recorded
source flux, rather than independently recomputing the whole nonlinear source.
The separate continuous ODE above checks that coupled trajectory in a fixed
tract. A complete energy identity for the moving planar graph is still not
claimed.

The closing refinement remains about 51.96 Pa then 39.54 Pa. The added internal
throat drag has not resolved that sensitivity. Losses elsewhere and pressure
trajectory accuracy remain separate questions; no damping constant has been
tuned to make this check look converged.

The twelve runs contain 7,171 accepted steps. Independent geometry, graph and
matrix-exponential checks cover 3,079 geometry snapshots at the first refinement.
Maximum discrepancies are 1.23e-8 Pa in graph pressure, 9.21e-15 Pa s in viscous
memory momentum and 1.28e-17 J in step viscous dissipation. Combined stored-air
residual is at most 7.07e-20 m3 per interval. The finest two closing intervals
dissipate 4.53e-7 J through these viscous segments; the two opening intervals
dissipate 7.16e-7 J. Their throat segment metrics are fixed, so their viscous
geometry work is zero. Changing-metric work is checked separately by the kernel.

## Location of the remaining refinement error

The [endpoint diagnostic](viscous_refinement_location.json) compares the 512 and
1,024 subdivision results at identical final geometry. The largest pressure
difference is in node 5, the last lip storage section, in both controls:

| Region | Closing maximum difference | Opening maximum difference |
|---|---:|---:|
| Mouth nodes | 0.00447 Pa | 0.03356 Pa |
| Throat nodes | 0.11337 Pa | 0.15412 Pa |
| Last lip node | 39.53881 Pa | 2.75778 Pa |

Projecting the charge and branch-flow difference onto the frozen inviscid
tract's wave modes localizes it further. In closing, an 81.30 kHz mode carries
99.89% of that projected difference energy; 95.09% of its pressure-energy
coordinate is at the last lip node. Its pressure contribution there is
39.53136 Pa. In opening, a 62.93 kHz mode carries 70.30% of the projected
difference energy and contributes 2.69916 Pa at that node.

This basis omits the source and loss-memory feedback and excludes circulation
invisible to node divergence. It is a diagnostic of the tract discrepancy, not
the complete damped-system spectrum or evidence that the error is inaudible.
It gives the next investigation a concrete location: the small, moving lip-end
storage and its presently missing local wall losses. Further throat damping is
not supported as a remedy by this comparison.

Seventeen focused CTests pass in Release and seventeen under ASAN/UBSAN. All
four actual moving viscous/thermal/radiation intervals pass under ASAN/UBSAN at
256 subdivisions/ms. Nineteen historical adapter outputs reproduce their full
JSON exactly, including the preceding thermal path.

## Reproduction

```sh
cmake --build build --target ducklips_slit_viscous_probe \
    ducklips_graph_continuous_glottis_probe ducklips_planar_inertive_motion_probe
python tools/slit_viscous_reference.py
python tools/continuous_viscous_reference.py
python tools/planar_continuous_reference.py --radiation \
    --thermal-slit-width 0.03 --viscous-slit-width 0.03
python tools/radiation_regressions.py --output out/viscous-regressions \
    --include-radiation --include-thermal
python tools/viscous_refinement_location.py
```

Use the existing scientific Python environment with NumPy, SciPy, mpmath and
Shapely. The [validation receipt](continuous_viscous_validation.json) pins the
fresh checks and their artifacts. No speech recording or intelligibility score
has changed. Performance work remains parked.

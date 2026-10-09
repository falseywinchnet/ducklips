# Glottal passage inertia and losses from the moving chambers

`GraphGlottalStorage::prepare_articulated` replaces the three fixed passage
coefficients with inertia derived from stored air and viscous resistance
derived from the same chamber dimensions. It also includes the kinetic-energy
contribution to mass-flow effort. The earlier fixed-coefficient path remains
a control. No new speech recording is attributed to this change.

The new component produces repeated fold opening/narrowing under a constant
700 Pa reservoir in a 100 ms control. No periodic tissue force, rest-opening
animation or supplied airflow drives that control. This is evidence of
emergent motion in these equations, not a validated human voice or a proven
stable physiological limit cycle. Its narrow openings and strong compression
need particular scrutiny.

## Derivation and choice of state

Each chamber has longitudinal thickness `h`, glottal width `w`, actual volume
`V`, aperture `A=V/h`, and ambient-reference mass volume `M`. Relative density
is `r=M/V`. A flow `Q` in reference-volume units has local velocity `u=Q h/M`.
The kinetic energy of half a chamber is therefore

```
K_half = (rho0 M/4) u² = rho0 h² Q²/(4 M)
L_half = rho0 h²/(2 M) = rho0 h/(2 r A)
```

At ambient density this becomes the familiar half-length air inertance.
At fixed mass, merely squeezing the chamber does not change its longitudinal
kinetic energy at fixed velocity: density and area change together. Updating
inertance as `rho0 h/(2 A)` would incorrectly assume ambient density throughout
that compression.

The three passages occupy the lower half-cell, the adjacent lower/upper
halves, and the upper half-cell. The middle passage has a common mass flow
through its two halves:

```
L0 = rho0 h_lower²/(2 M_lower)
L1 = L0 + L2
L2 = rho0 h_upper²/(2 M_upper)
G_j = 1/L_j
K = sum_j lambda_j² G_j(M)/2
Q_j = dK/dlambda_j = G_j lambda_j
```

`lambda` has units Pa s; it is the branch's potential/momentum state, not
literal total air momentum in kg m/s. In this half-cell representation it is
`rho0` times the length-weighted velocity integral along the passage.

Because `G` depends on the stored masses, moving air also carries
`kappa_i=dK/dM_i`. This is added to gas enthalpy in the passage momentum law.
For uniform area and equal flow in all passages, both chamber contributions
reduce to the Bernoulli kinetic head `rho0 u²/2`; the native probe checks that
limit and the integrated kinetic energy.

For fixed longitudinal thickness the kinetic energy has no explicit
volume derivative at fixed `M,lambda`. There is consequently no extra lateral
kinetic force to add to the fold pressure in this representation. Fold
pressure still derives from gas storage; kinetic head changes the solved
mass transfer and therefore pressure. This distinction is why changing
inertance alone would omit work, while adding an arbitrary closing force
would double-count a different mechanism.

The use of kinetic energy and relative enthalpy in compressible vocal-fold
coupling has precedent in [Mora et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8210679/).
The half-cell series discretization and finite-step implementation here are
our own restricted model, not a reproduction of that paper's equations.

## Finite-step work and losses

The symmetric discrete kinetic efforts satisfy

```
Delta K = sum_j Q_bar_j Delta lambda_j + sum_i kappa_bar_i Delta M_i
Q_bar_j = (G_old_j+G_new_j)(lambda_old_j+lambda_new_j)/4
```

`kappa_bar` uses symmetric coordinate secants of the inverse inertance and
averages the squared endpoint momenta. The middle inverse inertance is a
rational function of the two masses; its secants are evaluated analytically
without subtracting nearly equal values. These efforts enter the same coupled
three-flow solve as the gas/tissue and tract load. Passage momentum is not
reset when its metric changes. Combined source/tract work retains the balance
in [the previous coupling report](GRAPH_GLOTTAL_STORAGE.md).

The dissipative model uses the fully developed parallel-plate laminar law;
see the [MIT Couette/Poiseuille reference](https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_couet_and_pois/).
For a half-cell, conversion from actual volume flow to reference-mass flow
requires both the flow and its conjugate effort to change:

```
R_half = 6 viscosity h w²/(A³ r²)
R_passages = [R_lower, R_lower+R_upper, R_upper]
```

At constant density, doubling aperture reduces resistance eightfold. At
fixed mass, halving volume doubles resistance. Native controls check both
limits. The default terminal control has loss coefficient 1 multiplying
`rho0 h_upper² Q |Q|/(2 M_upper²)`. It is symmetric under flow reversal and
stands for lost exit kinetic head. It is not a jet separation/reattachment
model, and no claim of fitted anatomical accuracy is attached to it.

## Evidence

[Short-control receipt](glottal_geometry_reference.json):

- 241 kinetic-energy cases checked against 100-digit arithmetic, including
  nearly unchanged masses and signed momenta. Two invalid-mass cases reject.
  Maximum scaled chain-rule error is `9.99e-16`; maximum relative error in
  discrete mass effort is `2.84e-16`.
- Four independent continuous controls (forward, reverse, undriven decay,
  changing rest opening), each with eight native refinements over 5 ms.
  All 32 trajectories pass. Finest error-reduction factors are 3.866–4.005;
  the worst endpoint tract-pressure error is `2.86e-4 Pa`.
- The largest accumulated work residual across those trajectories is
  `1.37e-15 J`; largest stored-mass residual is `6.04e-21 m³`.
- The existing 32 fixed-passage continuous comparisons were rerun after the
  coupled solve's relative tolerance was tightened. Four articulated native
  controls and the kinetic-law probe pass ASan/UBSan. Seven focused CTests
  pass. Native stepping remains allocation-free.

The first 100 ms native control exceeded the probe's absolute passage-law
residual limit (`1.21e-5 Pa` versus `1e-5 Pa`). Tightening the coupled Newton
relative tolerance from `1e-10` to `1e-12`, retaining the `1e-8 Pa` absolute
term, removed that solver error. The physical parameters and the probe's
acceptance limit were unchanged.

The longer comparison also showed that small work residuals did not mean
accurate timing. At 19,200 steps over 100 ms, endpoint tract-pressure error
was about 167 Pa. At 153,600 steps it was still 8.78 Pa. These coarse results
are retained in the [long-control receipt](glottal_geometry_long_reference.json).
At 2,457,600 steps the doubling audit demonstrates second order (final
error reduction 3.998), but its normalized endpoint error `1.076e-4` still
misses the unchanged `1e-4` accuracy criterion. That failed flag is retained.
An [additional endpoint check](glottal_geometry_endpoint_reference.json) at
3,000,000 steps passes with normalized error `7.22e-5`, pressure error
`0.0229 Pa` and position error `2.36e-10 m`. Its observed order relative to
the previous point is 1.9997. It reuses the already step-refined continuous
oracle after verifying the source/binary hashes; it is not another ODE run.

A 33.3 ns step is a reference-accuracy result, not a practical rendering
representation. The work balance alone would have concealed how much
coarser trajectories differ in this near-contact regime.

## Physical limits exposed by the motion

A separate 76,800-step motion probe records sampled chamber areas and flows,
while collecting extrema at every step. It finds a minimum area of
`2.19e-8 m²`, a maximum local Mach number of 0.237, and density ratios from
0.729 to 2.212. Hydraulic Reynolds number reaches approximately 3,851.
These numbers belong to that numerical control; they do not establish that
its simple laminar loss law remains a valid human-glottis approximation.

The positive-volume model has no contact transport and no tissue-collision
branch. Very small remaining air volumes can develop large compression
forces. This requires resolving gas escape and tissue contact rather than
hiding the pressure with added damping or a volume floor. The model also
omits transverse fluid kinetic energy, separated jets, vortex shedding and
unsteady boundary-layer losses. The tract load is still the linear acoustic
graph, and moving mouth/lip replay is not part of these controls.

The sampled motion shows repeated openings under constant drive, but it is
not an acoustic realism score or a new NATO listening result. Performance
optimization remains parked; the fine-step reference is not being offered
as a production speech renderer.

## Reproduction

```sh
cmake --build build --target ducklips_graph_glottal_probe ducklips_glottal_kinetics_probe ducklips_glottal_geometry_motion_probe
python tools/glottal_geometry_reference.py
python tools/graph_glottal_reference.py
python tools/glottal_geometry_reference.py --duration .1 \
  --steps 19200 38400 76800 153600 307200 614400 1228800 2457600 \
  --modes forward --output out/glottal-geometry-long \
  --report reports/glottal_geometry_long_reference.json
build/ducklips_graph_glottal_probe --trajectory 3000000 0.1 forward --articulated
build/ducklips_glottal_geometry_motion_probe 76800 out/glottal-geometry/motion-76800.csv
```

The long audit records failed acceptance flags as well as its data before
returning failure. Source, binary and output hashes accompany the receipts;
[supplemental validation](glottal_geometry_validation.json) records sanitizer,
CTest and motion-probe provenance.

# Contour sealing and reopening with viscous relaxation

Follow-up: [the shrinking-pocket diagnosis](CONTOUR_POCKET_DIAGNOSIS.md)
examines the sustained-drive stop below with independent continuous equations
and isolates the role of the passage density/loss representation.

The contour source now completes two free-motion contact controls: the lower
cap seals and reopens in one, and the upper cap does so in the other. Tissue
motion is solved from initial velocity, elasticity, gas forces and damping;
no trajectory or event time is prescribed. Both chamber air masses remain
represented throughout. A closed cap has exactly zero flow and momentum.

This does not establish sustained voicing. A 700 Pa pressure-driven run
passes the first lower-cap seal, then stops near 10.91 ms as the lower
chamber becomes extremely small and compressed. The intended 100 ms run
is incomplete. There is no new speech recording, pronunciation result or
performance claim. The latest human feedback still belongs to the spatial-
pressure recording: November improved; Foxtrot remains jammed around /tr/.

## Why the time update needed to change

For fixed inverse inertia G and viscous resistance R, momentum obeys
`d lambda / dt = -R G lambda`. Let `x = dt R G`. The symmetric update has
multiplier `(1-x/2)/(1+x/2)`, which approaches -1 as x grows: a rapidly
relaxing physical mode becomes a persistent numerical sign alternation.
The cap's relaxation time shrinks toward zero near sealing.

A full backward update suppresses that stiff mode and permits contact, but
its first-order timing error is too large in three of the four tested open
controls at the finest tested grid. That failed accuracy result is retained
in [the backward comparison](contour_transition_reference.json), implemented
at commit `6b3a16d`. It is not the accepted accuracy result.

The new `ContourUpdate::relaxation` uses a continuous weight derived from
that relaxation equation, separately for each passage:

```
x    = dt R_endpoint G_endpoint
beta = coth(x/2) - 2/x
```

The weight tends to zero with x and to one for stiff relaxation. A series
near zero avoids cancellation. In the fixed-metric linear viscous case,
the resulting multiplier is exactly `exp(-x)` up to floating-point error.
There is no fitted switching threshold or new body damping parameter.
For moving geometry, midpoint and endpoint loss coefficients are blended
with the same weight. Exact scalar decay does not imply exact nonlinear
coupled dynamics; those are checked separately under time refinement.

## One finite work identity

For each passage, with endpoint inverse inertias G0,G1 and momenta l0,l1:

```
Qbar = (1-beta)(G0+G1)(l0+l1)/4 + beta G1 l1
C    = (1-beta)(l0*l0+l1*l1)/4 + beta l0*l0/2
N    = beta G1 (l1-l0)^2/2 >= 0
Delta K = Qbar Delta l + C Delta G - N
```

The finite mass and shape forces distribute `C Delta G` over the same
forward/reverse coordinate paths as the existing contour kernel. The gas
volume work and elastic contact work retain their exact finite identities.
All these terms enter the coupled tissue, stored-air and tract solve.

The diagnostic field `numerical_loss` records N. In the fixed-metric linear
viscous case, N supplies the within-step flow-variance part of viscous loss
omitted by `dt R Qbar^2`; together they reproduce the exact scalar energy
loss. Thus N is not extra artificial damping in that special case. For
variable geometry and nonlinear losses it is a nonnegative discrete balance
term, not a general estimate of physical dissipation or discretization error.
The `physical_loss_J` diagnostic records the explicit passage-loss and tissue-
damping terms; the report keeps N separate.

A sealing endpoint imposes G1=0, l1=0 and Qbar=0 for that cap. Its old kinetic
energy participates in the geometry work; it is not erased outside the ledger.
The inactive passage supports a constraint effort with zero flow. Reopening
introduces momentum from the coupled flow solve. Both adjacent gas states
remain present. There are no area floors or chamber-mass resets.

These are finite-step contact transitions, not an event-time root finder.
Recorded first closed/reopened times are sampled endpoints with resolution
one time step. Vanishing *chamber volume* is a separate unsupported condition.
The symmetric fixed-topology mode remains the default reference API behavior;
the new mode is selected explicitly by the research controls.

## Evidence and its limits

[The relaxation reference](contour_relaxation_reference.json) records:

- 160 weighted finite transfers against independent 80-digit calculations.
  Maximum scaled force discrepancy is 2.42e-13; the energy-chain discrepancy
  is 2.05e-15. Cases include signed contours and sealing/reopening pairs.
- 98 relaxation-rate fixtures from zero through 1e12. Maximum weight error
  is 4.62e-16; fixed-passage decay error is 1.84e-16 absolute.
- Four fresh coupled open-motion refinements against the preserved independent
  continuous endpoints from commit `5bc0d1ace9273b158346921fd1a652edbb8c508f`.
  At 15,360 steps over 5 ms, scaled endpoint errors are below 1.00e-6;
  final refinement ratios range from 3.75 to 3.99. The continuous oracles
  are reused and their source/definition bytes checked, not reintegrated.
- Two free-motion contact cycles on six grids, from 500 to 16,000 steps over
  2 ms. Each contains exactly one sealing and one reopening. Finest successive
  scaled endpoint differences are 1.82e-5 (lower) and 1.42e-6 (upper).
  These are self-refinement checks, not independent full-cycle ODE solutions.
  Contact alignment affects ratios; no universal second-order cycle claim
  follows from these two controls.

At the finest grid, the lower cap's first closed/reopened endpoints are
77.875/870.375 microseconds, and the upper cap's are 86.625/520.125 microseconds.
Maximum represented cap energy immediately before sealing is below 1.43e-28 J.
The accumulated N terms are 5.80e-13 and 5.41e-13 J, respectively, and decrease
under the final refinements. Maximum step work residual is below 1.36e-18 J;
maximum step mass residual is below 3.84e-23 cubic metres of reference-density
air. These conservation checks supplement rather than replace trajectory checks.

Twenty focused CTests pass. Fresh native symmetric open controls still pass.
Address/undefined-behavior sanitizer runs pass the weighted kernel, coupled
controls and both contact cycles; weighted-transfer and relaxation-weight
fixture output is byte-identical to native. [Validation receipts](contour_relaxation_validation.json)
record source/binary/output hashes, commands, logs and the incomplete drive.

## Remaining pressure-driven failure

For 76,800 planned steps over 100 ms, the solve stops after 8,376 accepted
steps (10.90625 ms). With four times as many planned steps it stops after
33,505 steps (10.9065755 ms). Both runs have one lower-cap seal and no reopening.
At the finer stop, lower chamber volume is 8.04e-16 cubic metres, relative
density 2.59, and pressure about 279 kPa above ambient. The upper chamber
remains about 1.37e-10 cubic metres. The time step was 0.326 microseconds.

These stops do not establish whether a resolved continuous solution bounces
under compression, empties the chosen control volume, or needs a different
control-volume partition. That distinction needs investigation. Removing a
cell would require explicit mass, momentum and energy transport. The current
solver rejects loss of positive chamber storage; it does not silently delete
the retained air. Sustained voicing, connection to the speech renderer, the
Foxtrot /tr/ correction and realistic NATO speech remain unfinished.
Performance optimization remains parked.

Reproduce with the scientific Python environment used by the repository:

```
cmake -S . -B build
cmake --build build
python tools/contour_relaxation_reference.py
python tools/contour_open_regression.py
build/ducklips_contour_contact_probe 76800 --relaxation
build/ducklips_contour_contact_probe 307200 --relaxation
```

The contact diagnostic intentionally exits successfully when it records a
solver stop; inspect its JSON `status` and accepted/planned step counts.

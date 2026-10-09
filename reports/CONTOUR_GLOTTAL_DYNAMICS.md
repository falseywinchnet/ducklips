# Coupled open-contour glottal dynamics

The medial contour now drives a working gas/tissue/passage/tract solve.
`ContourKinetics` derives both mass and shape forces from the same kinetic
energy. `ContourGlottalStorage` solves the four tissue displacements and
three passage flows together; `GraphContourGlottal` connects its outlet to
the existing total-flow acoustic graph.

This report records the initial open dynamic reference. The subsequent
[sealed-cap extension](SEALED_GLOTTAL_DYNAMICS.md) adds elastic contact and
fixed-topology sealed phases; automatic sealing/reopening remains unfinished.
This is a dynamic reference component, not a new speech recording.
The latest human pronunciation feedback still belongs to the spatial-pressure
recording: November improved; Foxtrot remains jammed around /tr/. Performance
optimization remains parked. These correctness controls carry no realtime claim.

## What changed physically

The previous rectangular model made a whole chamber's volume proportional to
one narrow control opening. The contour instead integrates the air beside
both tissue controls. Its chamber volumes depend on both openings, and its
pressure locations follow their air-volume centroids. This changes both gas
pressure work and the inertia of the air paths.

The profile, widths, tissue masses, springs, damping and ideal-gas law retain
the assumptions documented in [the geometry reference](GLOTTAL_CONTACT_GEOMETRY.md)
and [the storage reference](GLOTTAL_STORAGE.md). No damping adjustment,
aperture floor, pressure clamp, or prerecorded excitation was added.

Let `M_i` be chamber air mass divided by ambient density, `V_i(a,b)` its
volume, and `r_i=M_i/V_i`. Each of the three paths is split at the fixed
chamber boundary. For its portions in each chamber:

```
I1_ji = integral dz/A(z)
I3_ji = integral dz/A(z)^3
L_j   = rho0 sum_i [(V_i/M_i) I1_ji]
G_j   = 1/L_j
K     = (1/2) sum_j lambda_j^2 G_j
Q_j   = G_j lambda_j
R_j   = 12 viscosity width^2 sum_i [(V_i/M_i)^2 I3_ji]
```

`Q` is reference-density mass volume per second, not local physical volume
flow. The optional outlet kinetic-head loss retains the earlier symmetric
terminal-loss assumption. It is not a jet-separation law.

At fixed tissue coordinates, the mass effort is `dH/dM_i + dK/dM_i`.
At fixed masses and canonical passage momenta, the area effort is

```
E_area,k = dK/darea_k - sum_i P_i dV_i/darea_k.
```

Each left/right tissue mass at that control gets fluid force
`-width * E_area,k`. The rest-area actuator contributes
`sum_k E_area,k * d(rest_area)`. Thus moving centroids, changing density,
changing inlet/outlet geometry and gas displacement all participate in one
energy account. Replacing only chamber volume would omit the kinetic shape
force and its actuator work.

This is still a reduced model: density is uniform within each chamber,
passage velocity follows its local area, and shear loss is the existing
quasisteady parallel-plate assumption. Energy consistency establishes the
accounting of these assumptions, not their physiological adequacy.

## Finite motion without subtracting nearly equal energies

The kinetic update uses a symmetric coordinate-increment gradient over
`(M_lower,M_upper,a,b)`, averaged over both momentum endpoints. Its passage
flow is

```
Qbar_j = (G_old,j + G_new,j) (lambda_old,j + lambda_new,j) / 4.
```

It satisfies the finite chain rule

```
Delta K = sum_j Qbar_j Delta lambda_j + sum_k Ebar_k Delta coordinate_k.
```

The internal paired scalar carries endpoint values and their divided
difference through products, quotients, logarithms and clipped path endpoints.
It becomes an ordinary derivative at coincident endpoints. A convergent
series handles nearly uniform area; well-separated coordinate endpoints use
ordinary secants to avoid cancellation in long composed expressions. There
is no numerical perturbation used to define the constitutive force.

Open-contour chamber volumes are affine in the two areas, so their existing
volume gradients are also exact finite volume differences. The gas uses its
previously checked mass/volume discrete gradient. Tissue springs and inertia
use midpoint work; positive viscous, terminal and tissue damping terms enter
as dissipation. The coupled Newton solve uses a numerical Jacobian of these
residuals, bounded iterations and line search. It works with tissue
**displacements**, avoiding loss of precision from subtracting positions.

The acoustic graph is eliminated into its free average load plus its
unit-flow response, then restored with the accepted outlet flow. The source
ledger is

```
Delta E_source + tissue_damping + passage_loss + work_to_tract
  - reservoir_work - rest_actuator_work - external_tissue_work = 0.
```

Failure preserves public state and diagnostics. Stepping allocates no heap
storage. The old rectangular source remains available as a separate control.

## Independent checks

[The initial reference report](contour_glottal_reference.json) contains:

- 180 arbitrary-precision instantaneous energy/force/resistance checks and
  60 finite-transfer checks at 80 decimal digits. Largest scaled force/transfer
  discrepancy: `1.84e-14`.
- 1,000 native randomized finite-chain, coincident-limit and time-reversal
  checks in the CTest kernel; largest scaled work discrepancy `1.02e-14`.
- Four 5 ms continuous controls: forward flow, reverse flow, undriven pressure
  relaxation and moving rest opening. Their independent DOP853 equations use
  complex-step derivatives of separate interval-integration code, not the
  native finite-change forces. Refining its maximum step from 2 to 1 microsecond
  changes normalized endpoints by at most `3.40e-10`.
- Six native time grids per control, 120 through 3,840 steps. Three controls
  pass the stipulated second-order convergence and endpoint-accuracy checks.
  Relaxation meets endpoint accuracy but **fails the coarse-grid ratio check**;
  this failure remains in that report.

[The additional endpoint report](contour_glottal_endpoint_reference.json)
reuses the verified continuous relaxation endpoint and checks 7,680 and
15,360 native steps. Error decreases by factors **3.8776 and 3.9222**. Its
finest normalized error is `5.53e-7`, and maximum tract-pressure endpoint
error is `0.000312 Pa`. All four controls now meet the original criteria.
This is additional native refinement, not another continuous-oracle run.

A separate asymmetric-body control drives unequal external tissue forces and
moving rest area against a nonzero outlet effort/impedance. Its accumulated
work discrepancy is `2.00e-16 J`; maximum force residual is `1.24e-14 N`.
Zero-drive equilibrium, rejected input, state restoration and no-allocation
checks pass. Eleven focused CTests pass, including adjacent older components.
Address/undefined-behavior sanitizer runs pass for the new kernels, four
coupled controls and asymmetric control. Both sanitized arbitrary-precision
fixture batches are byte-identical to native outputs.

## The next boundary: contact

Under the same standard body and 700 Pa reservoir, the open-contour solver
stops as its lower opening approaches zero. These are the **last accepted
states before failure**, not certified contact events and not completed
100 ms trajectories:

| Planned steps over 100 ms | Last accepted time, ms | Lower opening, m² | Lower volume, m³ | Maximum density/ambient before stopping |
| ---: | ---: | ---: | ---: | ---: |
| 19,200 | 10.604167 | 1.506e-8 | 1.413e-9 | 1.00548 |
| 76,800 | 10.605469 | 5.568e-9 | 1.389e-9 | 1.00570 |
| 307,200 | 10.606120 | 1.105e-9 | 1.377e-9 | 1.00572 |

At the finest last state, lower area is shrinking at `0.00679 m²/s`;
linear extrapolation would cross zero inside the next time step. Both chambers
retain positive volume. Their density ratios are `0.99810` and `0.99839`;
relative gas pressures are about `-265` and `-225 Pa`. This demonstrates a
trajectory approaching a sealed cap with retained wedge storage. It does
**not** show that the earlier 100 ms rectangular compression problem has
been cured over a complete cycle: this new solver has not crossed contact.

The geometry component can describe contact, but this dynamic solver
explicitly requires both control openings positive. Tissue collision,
passage momentum at sealing/reopening, and disappearing-cell state transport
still need their own consistent dynamics. The short trajectories verify the
open equations. They do not establish perceptual realism, physiological
calibration, utterance integration or production efficiency.

## Reproduce

Use the project's existing scientific Python environment
`/Users/joshuahkuttenkuler/miniforge3/bin/python`:

```
cmake -S . -B build
cmake --build build --target ducklips_contour_kinetics_probe ducklips_contour_glottal_probe ducklips_contour_source_controls ducklips_contour_contact_probe
python tools/contour_glottal_reference.py
# The initial relaxation-ratio assertion fails, intentionally retained above.
python tools/contour_glottal_endpoint_reference.py
build/ducklips_contour_source_controls
build/ducklips_contour_contact_probe 19200
build/ducklips_contour_contact_probe 76800
build/ducklips_contour_contact_probe 307200
```

The endpoint script verifies the original report's source/binary hashes before
reusing its oracle. [Validation receipts](contour_glottal_validation.json)
record the contact outputs, focused test log and sanitizer evidence separately.

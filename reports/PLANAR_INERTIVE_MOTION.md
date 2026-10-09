# Air momentum in the moving tongue-contact geometry

The moving polygonal volumes and openings now drive the existing canonical
pressure-and-momentum integrator used by the current speech research path.
This is an inviscid diagnostic on the five-region partition, not the old modal
renderer and not a new speech recording. Closure and release remain two
separately initialized 1 ms intervals.

The preceding quasi-steady control's 12% density reduction was strongly
affected by its arbitrary resistance. Replacing that control with explicit
air inertia reduces the worst sampled reduction to 3.7%. The remaining front
chamber pressure drop is consistent with its expanding volume and deliberately
sealed exterior. Neither result establishes realistic human pressure or
resolves the audible Foxtrot /tr/ issue.

## What the 12% result was measuring

During the saved 5.062–5.063 s release interval, the small lower-front-left
region grows from 2.2390 to 2.5495 mL: 0.3105 mL, or 13.9%, in 1 ms. Its
opening into the neighbouring front region is broad, growing from roughly
439 to 486 mm². The old fixture's chosen Darcy mobility nevertheless admits
only 0.00679 mL net reference-volume air into that region over the interval.
It behaves almost as if separately sealed.

The new control admits 0.21585 mL into it. Its final pressure changes from
−16634.3 Pa in the resistance control to −5182.8 Pa with momentum. The two
other front regions finish at −4646.6 and −4681.9 Pa. At exactly the same total
final front-chamber mass and volume, uniform pressure would be −4747.7 Pa.
Thus the residual drop has a substantial whole-chamber component, rather
than being only a numerical isolation of the smallest region.

All three front regions together expand by 0.5280 mL while the exterior is
sealed. The developing tongue opening admits about 0.03719 mL net from the
rear. These are reference-volume mass units under the linear storage law.
They are not measurements of a person's mouth. Opening a physical lip port
will change the boundary-value problem; adjusting an internal resistance
cannot substitute for that missing boundary.

## Coupled storage and momentum

Let `B` denote the oriented region/branch incidence matrix, `q = C p` the
excess stored air, and `lambda` the branch momentum coordinate. The diagnostic
advances the same explicitly assumed canonical equations as `MovingTreeReference`:

```
C = V / (rho c²)
G = A / (rho d)
Q = G lambda
q' = -B Q - V'
lambda' = Bᵀ p
```

Here `A` is the actual open grid-face area, `d` the normal grid-centre distance,
and `rho = 1.14 kg/m³`, `c = 350 m/s`. On a straight uniform duct this is the
usual lumped air inertance `rho d/A`. There is no viscous resistance assigned
to an artificial grid face. A coarse irregular cut region does not establish
a uniform flow path of that length; spatial convergence and better geometric
inertance are still required there.

Each two-stage Gauss step builds its stage volumes and openings directly from
the moving native contours. It does not interpolate endpoint metrics. The
wall-volume source is held at `-(V_end - V_start)/dt`, making endpoint stored
mass `V + q` conservative. This is an approximation to the within-step wall
velocity, so fourth-order homogeneous acoustic integration does not imply
fourth-order accuracy for general prescribed wall motion.

Actual opening/closure times split the time steps. Zero-inertive-conductance
branches carry zero flow, and their latent momentum is cleared as in the
canonical speech closure handling. No pressure averaging occurs on reopening.
The tool's union graph is a tree; this does not add arbitrary-cycle acoustic
support or solve general moving region identities.

`MovingTreeReference::average_flow` now exposes the actual Gauss-weighted
branch flow. Its product with the step duration is the transported volume;
it is also conjugate to a held series-pressure source. The core state update
formula is unchanged. A new source-work check verifies that transported volume
and supplied acoustic energy agree. The output vector is allocated during
preparation, not during the integrator's step. The offline geometry/reporting
tool itself is not an allocation-free renderer.

## Checks against independent calculations

Ten full moving runs use 8, 16, 32, 64 and 128 subdivisions per event, plus
the contact times: 506 pressure/momentum steps. Independent Boolean geometry
checks 1,528 endpoint and Gauss-stage geometries. A separate dense collocation
solve evolves all charges and momenta from initial conditions using independent
volumes and apertures. It never resets to intermediate native states.

| Maximum discrepancy over the moving runs | Result |
|---|---:|
| Region volume against Boolean geometry | 4.68e-19 m³ |
| Port area | 5.43e-18 m² |
| Endpoint pressure against independent solve | 2.67e-9 Pa |
| Branch momentum against independent solve | 3.16e-13 Pa s |
| Average transported flow against independent solve | 3.94e-15 m³/s |
| Local stored-air balance | 6.47e-21 m³ |
| Whole-domain stored-air balance | 1.36e-20 m³ |
| Sealed rear/front group balance | 6.78e-21 m³ |
| Flow through a fully closed throat in checked intervals | exactly zero |

The existing conservative continuation-witness check matches all five regions
over every step. Its quarter-time point checks supplement the swept-boundary
construction; they are not a complete detector of arbitrary transient topology.
The stored contact endpoint uses actual native sampled coordinates for the
independent Boolean audit. Stage points are away from contact and independently
reconstructed from the original corresponding vertices.

The 64-to-128 subdivision endpoint differences are 0.0326 Pa for closure and
0.1843 Pa for release. These compare successive discretizations, not an exact
continuous solution or a spatially resolved mouth.

Ten analytical controls cover two distinct physical limits:

- Nine closed uniform-duct runs use 8, 16 and 32 spatial cells at three time
  resolutions. The pressure follows the exact discrete standing wave, with
  roughly sixteenfold time-error reduction per doubling and less than
  1.45e-5 Pa error at the finest step. Relative energy drift is at most
  6.67e-16. The first-mode frequency approaches the continuum value
  1458.333 Hz: 1448.981, 1455.992 and 1457.748 Hz as the mesh is refined.
- A sealed expanding single region matches the exact declared storage law,
  `p = rho c² (M_initial/V - 1)`, within 1.69e-11 Pa. This verifies that wall
  motion changes pressure while preserving its stored air.

These controls test wave propagation, source work and moving storage. The
moving dense solve checks the discrete implementation; it does not separately
validate the chosen governing equations.

## Physical limits and next integration work

The five-region model has no viscous, thermal, jet or radiation losses. The
maximum sampled normal face velocity during release is 46.4 m/s (Mach 0.133)
under this model. A lossless linear pressure-and-momentum control cannot
validate that flow as a human constriction jet. The body/blade union is also
not constrained to constant total tissue volume, and the fixed 22 mm extrusion
is not individually measured anatomy.

Depth-aware planar acoustics has precedent, but a depth parameter alone does
not prove equivalence to a three-dimensional vocal tract; the
[2.5D vocal-tract study](https://arxiv.org/abs/1909.09585) explicitly develops and
evaluates a particular depth-aware wave model. This diagnostic has not adopted
or reproduced its acoustic validation. Likewise, measured vocal-tract losses
depend on constrictions and expansions; the
[Birkholz–Häsner study](https://www.isca-archive.org/interspeech_2024/birkholz24_interspeech.html)
does not justify drag at arbitrary numerical cell boundaries.

The next required boundary work is to retain actual glottal and lip openings
in the planar domain and couple them to source and radiation states. Physical
neck loss, nonlinear mean flow, spatially adequate inertance, continuous full
trajectory state transport, nasal routes and acoustic comparison remain
unfinished. The latest human feedback still belongs to the earlier
spatial-pressure recording; there is no new NATO score. Performance work is
still deferred.

## Reproduction and provenance

```
cmake --build build --target ducklips_planar_inertive_motion_probe ducklips_moving_tree_probe ducklips_canonical_closure_tests
python tools/planar_inertive_motion_reference.py
ctest --test-dir build --output-on-failure -R 'ducklips_moving_tree_reference|ducklips_canonical_closure_tests'
```

The [receipt](planar_inertive_motion_reference.json) records commands,
source/binary/input hashes, per-run artifact hashes, analytical controls and
time refinement. Full native states and independent dynamics are in
`out/planar-inertive-motion/`.

The two affected native core checks passed fresh in Release and under
AddressSanitizer/UndefinedBehaviorSanitizer. Sanitizers also passed all ten
analytical controls and both actual moving events at 16 subdivisions, with
their independent audits; see `out/planar-inertive-asan/receipt.json`.
The full historical test suite was not rerun for the added diagnostic flow
output. No benchmark or new speech render is counted in this validation.

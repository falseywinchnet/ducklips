# Finite narrow-gap loss for the continuous contour

`ContourPressureDissipation` now integrates a symmetric five-rate resistance
matrix over the continuous contour. Its local conductance comes from an
offline gas-kinetic reference and a compact native interpolant. The moving-tip
control has a finite full-domain loss of **1.29766786 µW**, without a minimum
gap or omitted endpoint interval. The inherited continuum no-slip law diverges
for that same moving-tip control.

This is a constitutive component of the new continuous source. It does not
advance momenta, drive the tissue, complete a phonation cycle, or change any
NATO recording. The user's latest listening result remains November improved
and Foxtrot crowded around /tr/ on the spatial-pressure recording. Performance
optimization remains parked.

## Constitutive reference and numerical representation

The offline reference solves the linearized, isothermal BGK plane-channel
problem with fully diffuse stationary walls. [Sharipov and Seleznev's review](https://srd.nist.gov/jpcrdreprint/1.556019.pdf)
supplies the integral equation, reduced-flow normalization, six comparison
values in Table 1 column f (page 679), and the large-delta asymptote in equation
3.88 (page 691). BGK is a model collision operator; the review distinguishes
its results from the full Boltzmann equation. These are benchmark conditions,
not measured properties of a human vocal-fold surface.

For gap h, density rho and viscosity eta, define thermal speed vth and
rarefaction parameter delta by

```
vth = sqrt(2 R T)
delta = h rho vth / (2 eta).
```

Our cell-average reference integrates the spatial exponential transport
kernel exactly within each cell and integrates positive molecular velocities
on a logarithmic coordinate. The symmetric Toeplitz system is

```
(I - K) V = K 1 / delta
G(delta) = mean(V).
```

The positive-speed Gaussian weights integrate to one half. No empirical
normalization is applied. The reference uses 512 through 32768 spatial cells
as needed and compares 16- and 24-point velocity rules per log interval.
The first attempt stopped at 2048 cells and failed its original 3e-4 spatial
change criterion at delta 68.13 and 100. [That failure is retained](slit_kinetic_initial_grid_failure.json).
The revised calculation refines those grids and uses a stricter 1e-5 change
criterion throughout; it does not relax the failed threshold.

[The refined reference](slit_kinetic_reference.json) contains 201 log-spaced
points from delta 1e-8 to 100. Native code stores 200 cubic intervals for
log(G), keeping G positive. Analytic tail values and slopes join the table
continuously with continuous first derivatives. The tails are

```
small delta: G ~ [-log(delta) + 3(1-EulerGamma)/2] / sqrt(pi)
large delta: G ~ delta/6 + 1.0162 + 1.0653/delta - 2.1354/delta².
```

The small-delta expression follows by expanding the single-flight integral
and retaining its logarithmic term. The large-delta expression is the cited
BGK asymptote. The former is checked at 1e-4, 1e-6 and 1e-8; the latter at
100, 200 and 500 against three spatial refinements and a second-order
Richardson estimate. These checks support the chosen joins, not a uniform
rigorous error bound over every positive delta. The native table is generated
only after its reference assertions pass. No kinetic solve runs in the engine.

## Loss from the same local continuity field

Let r be density relative to ambient, A the local airway area, and w the
five rates `[Qin, Qmiddle, Qout, lower_area_rate, upper_area_rate]` used by
`ContourPressureKinetics`. It reconstructs `Q(z)=F(z) dot w`, where Q is
mass flow expressed as ambient-reference volume per second.

The local BGK mobility gives average axial velocity

```
v = -h G(delta) / (rho vth) * dp/dz.
```

The gas-storage effort mu satisfies `dmu=dp/r`, and `Q=r A v`.
With `h=A/width`, substitution yields

```
Q = -[r A² G / (rho0 vth width)] * dmu/dz
R'(z) = rho0 vth width / [r A² G]
Rij = integral R'(z) Fi(z) Fj(z) dz
power = w^T R w.
```

The same flow reconstruction therefore supplies flow/flow, flow/wall and
wall/wall loss terms. Each quadrature contribution is a positive scalar
outer product, so the matrix is positive semidefinite on admissible rates.
Closed exterior flows remain constrained to zero and are rejected by the
power API. Zero viscosity returns zero resistance as a limiting case.

For consistency with the present isentropic storage state, the local thermal
speed is `sqrt(2/gamma) c0 r^((gamma-1)/2)`. Viscosity remains the configured
constant. This uses the isothermal conductance locally; it does not add a
heat equation or thermal creep to the engine.

At large delta, `G ~ delta/6` recovers precisely the familiar local
resistance `12 eta width²/(r² A³)`. At a moving contact tip, continuity gives
`Q/(r A) -> tip_speed`. Hence the new loss density tends to

```
rho0 r vth width * tip_speed² / G(delta),
```

which vanishes logarithmically as A approaches zero. Its axial integral is
finite. This is the consequence of the chosen constitutive limit, rather
than an area floor or a fitted contact damping coefficient.

## Integrating to the actual endpoint

Native integration uses an endpoint coordinate `distance=span*t³` on a
contact segment, with its exact Jacobian. The full interval is retained.
Eight- and sixteen-point Gaussian rules adapt against a shared whole-segment
error budget; the endpoint contribution becomes small in absolute terms.
A solely local relative criterion cannot accept a logarithmic endpoint.
Recursion is bounded at 16, with failure returned if accuracy is unresolved.

This check exposed an existing cancellation issue in restricted storage
integrals. Recomputing very small area and shape-function values from rounded
axial coordinates could defeat relative quadrature near contact. Storage now
interpolates both endpoint shape values and endpoint areas within each
integration interval. Contact clipping and local flow sampling share the same
contact origin. No density, volume or area floor was added. Full storage,
inverse-moment and flow checks were rerun after this change.

An independent axial Gauss-Kronrod integration uses a different quadrature
and no native endpoint transform. In the constant-density moving-tip case,
it also uses an analytic flux rather than the native flow reconstruction:

| Diagnostic omitted distance | Remaining integrated BGK loss |
| --- | --- |
| 1e-5 m | 0.991301993 µW |
| 1e-6 m | 1.245132648 µW |
| 1e-7 m | 1.293512530 µW |
| 1e-8 m | 1.297388746 µW |
| 1e-9 m | 1.297647971 µW |
| zero, complete domain | 1.297667860 µW |

The native full-domain result is 1.297667862 µW. These diagnostic omissions
are used only to show convergence; none enters the engine calculation.

## Verification and remaining scope

[Independent checks](contour_pressure_loss_reference.json) establish:

- Six published BGK benchmark comparisons: maximum relative discrepancy
  5.10e-5 against the printed values.
- Refined table: maximum final spatial change 9.12e-6 and velocity-quadrature
  change 5.52e-13.
- 600 withheld quarter/mid/three-quarter interpolation points: native
  conductance differs from the direct kinetic reference by at most 1.71e-6.
- Twelve open/contact loss matrices: maximum scaled discrepancy 3.60e-10
  against independent axial integration; its tolerance refinement changes
  the matrices by at most 1.31e-11. Nonnegative eigenvalues hold within
  floating-point roundoff.
- Fresh 80-digit storage evaluation of 72 cases: maximum scaled error
  1.43e-15 after the endpoint representation change.

[Validation receipts](contour_pressure_loss_validation.json) also record fresh
flow, inverse-moment and prescribed-geometry seam regressions, three focused
CTests, address/undefined-behavior sanitizer tests, and byte-identical release
and sanitized loss/conductance fixture outputs. Source, generated table,
binary and output hashes identify the tested implementation.

The physical approximation remains material: a local, steady, linearized
parallel-wall BGK mobility is not a solution for the actual moving wedge,
finite channel, humid diatomic gas, thermal gradients or unknown wall
accommodation. Its finite contact limit resolves the continuum-law divergence
within this model; it does not establish physiological validity or speech
realism. Compatible generalized momenta, tissue dynamics, topology changes,
source/tract coupling and listening tests remain necessary.

Reproduce with the repository scientific Python environment:

```sh
OPENBLAS_NUM_THREADS=1 python tools/slit_kinetic_reference.py
cmake -S . -B build
cmake --build build --target ducklips_contour_pressure_loss_probe ducklips_contour_pressure_flow_probe ducklips_contour_pressure_probe
OPENBLAS_NUM_THREADS=1 python tools/contour_pressure_loss_reference.py
# The validation script also expects matching targets in .build/asan,
# configured with -fsanitize=address,undefined -fno-omit-frame-pointer.
OPENBLAS_NUM_THREADS=1 python tools/contour_pressure_loss_validation.py
```

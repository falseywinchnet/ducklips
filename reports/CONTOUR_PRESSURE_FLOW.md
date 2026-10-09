# Airflow and inertia from continuous local mass conservation

`ContourPressureKinetics` now reconstructs axial airflow from the continuous
gas field and its moving contour, then integrates its kinetic energy. The
result is a coupled five-rate inertia matrix, including airflow caused by
wall motion. It does not contain a fixed internal chamber seam.

The same reconstruction exposes a separate limitation of the inherited
viscous law: a moving contact point gives logarithmically divergent no-slip
slit dissipation. The new component records the divergence coefficient; it
does not replace the divergence with a finite cutoff. Compatible momenta,
a valid narrow-gap loss treatment and time integration are still needed
before this can replace the speech source. There is no new audio result.

## Reconstructing the flux rather than assigning it

The storage component supplies continuous density r(z), area A(z), two
mass moments M, compliance C=dM/du, and geometry derivative D=dM/da. Let the
five rates be

```
w = [Qin, Qmiddle, Qout, lower_area_rate, upper_area_rate]
Mdot = [Qin-Qmiddle, Qmiddle-Qout]
udot = C^-1 (Mdot - D adot).
```

Qmiddle is the integrated transfer associated with the two overlapping
mass moments. It is not flow at an internal chamber boundary. For the
fixed shape functions phi, local continuity requires

```
d(Ar)/dt = r (phi dot adot) + A r^(2-gamma)/B (phi dot udot)
dQ/dz = -d(Ar)/dt
Q(z) = Qin - integral_0^z d(Ar)/dt dz.
```

An equivalent expression integrates backwards from Qout. The implementation
uses the nearer end, avoiding subtraction of large cumulative flows near
a closed endpoint. Restricted storage integrals provide the coefficients
without duplicating the gas constitutive law. Occupied regions contribute
zero air and zero flux; a closed exterior port requires its flow to be zero.

Integration against the shape functions recovers

```
(1/ell) integral_between_tissue_centres Q(z) dz = Qmiddle,
```

where ell is the distance between those centres. Thus both the local and
moment forms of continuity describe the same field.

## Fluid motion includes wall-motion inertia

Write the reconstructed flux as Q(z)=F(z) dot w. Axial kinetic energy is

```
K = (rho0/2) integral_air Q(z)^2/[r(z) A(z)] dz
  = (1/2) w^T L w
Lij = rho0 integral_air Fi(z) Fj(z)/[r(z) A(z)] dz.
```

This symmetric matrix includes flow/flow, flow/wall and wall/wall terms.
Simply changing the old three independent passage inertias would omit these
couplings. L is positive semidefinite on admissible rates. It need not be
strictly positive: for example, some instantaneous wall motions can be
accommodated by local density changes with no axial flow. Tissue inertia
must participate in the full mechanical system.

For a uniform tube at uniform density and constant throughflow, the native
component recovers `K=rho0 H Q²/(2 r A)`. Near a sealed contact endpoint,
the admissible local flow vanishes and the axial kinetic energy remains
integrable. A forbidden nonzero closed-port rate is rejected; its formal
column is not included as a hidden, zero-energy degree of freedom.

For a future variational integrator, one useful coordinate choice is
`[J,M0,M1,a,b]`, with

```
Qin = Jdot
Qmiddle = Jdot - M0dot
Qout = Jdot - M0dot - M1dot.
```

This is a linear rate transformation of L. The full equations must account
for its dependence on mass and geometry and include tissue inertia and gas
energy. The present API computes the instantaneous metric; it does not yet
advance those generalized momenta or transport them through topology changes.

## Why the old viscous law cannot be carried over unchanged

Let a contact tip be zc and locally let A(z)=s(z-zc), with finite positive
contact density. Local continuity gives

```
Q/(r A) -> zc_dot
zc_dot = -A_t(zc)/A_z(zc).
```

The existing no-slip parallel-plate dissipation density would therefore obey

```
12 eta width² Q²/(r² A³)
    ~ [12 eta width² zc_dot²/abs(s)] / distance_from_contact.
```

Its integral diverges logarithmically when the contact point moves. The
density cancels from the leading coefficient. If the contact point is
stationary, this leading obstruction vanishes; this statement does not by
itself evaluate all remaining loss terms.

The native `contact_loss` matrix returns the quadratic coefficient in the
two area rates. It is explicitly a divergence diagnostic, not a finite
Rayleigh damping matrix or a newly introduced loss parameter. At a zero-area
plateau endpoint it follows the closed-side contour convention.

An independent moving-tip control has tip speed 0.1125 m/s. Integrating the
old loss only down to successively smaller diagnostic distances gives:

| Axial distance from tip omitted | Truncated viscous power |
| --- | --- |
| 1e-5 m | 1.03777e-6 W |
| 1e-6 m | 1.51585e-6 W |
| 1e-7 m | 1.99393e-6 W |
| 1e-8 m | 2.47202e-6 W |
| 1e-9 m | 2.95010e-6 W |

Every decade adds approximately 4.78082e-7 W. The measured coefficient per
natural logarithm is 2.076283125e-7 W, matching the derived and native
coefficients. No such cutoff enters the engine component.

A gas-wall treatment with an appropriate small-gap limit needs investigation.
Maxwell slip is a possible reference in its valid regime, not an automatic
justification for extending a continuum law to zero gap. [Lockerby's
first-order slip analysis](https://www.cambridge.org/core/journals/flow/article/firstorder-analysis-of-slip-flow-for-micro-and-nanoscale-applications/490BD4648720180653EFF9F1733CD513)
relates slip length to mean free path and accommodation, and explicitly
restricts the continuum slip treatment to low rarefaction. It identifies
additional gas-kinetic effects at higher rarefaction. No slip law or fitted
minimum gap is adopted in this change.

## Verification and scope

[The independent flow audit](contour_pressure_flow_reference.json) covers
12 open, lower-contact and upper-contact geometries, including the old seam
and exact zero-cap endpoints. Separate SciPy Gauss-Kronrod integration
reconstructs the field and inertia. Maximum scaled discrepancies are
1.48e-15 for the inertia matrix and 5.51e-14 for local flow coefficients.
Integrated moment continuity is satisfied within 6.67e-16. The tested inertia
matrices are positive semidefinite within floating-point roundoff. Two
reference tolerances produced identical integrals in these cases; this is
not a time-integration convergence claim.

The storage refactor exposes restricted intervals but preserves all full-
storage fixture outputs, inverse-moment results and the seam sweep byte for
byte against commit `980c679`. Those original independent 80-digit oracles
are reused, not regenerated. Native uniform-flow/energy controls, both
contact supports and closed-flow rejection pass. Two focused CTests and
address/undefined-behavior sanitizer controls pass; sanitized flow fixture
outputs are byte-identical to release outputs.

[Validation receipts](contour_pressure_flow_validation.json) retain commands,
source/binary/input/output hashes, regression provenance and sanitizer logs.
There is no physiological calibration, sustained-voicing result, new NATO
listening result or production-speed claim in these component checks.
Performance optimization remains parked.

Reproduce with the repository's scientific Python environment:

```
cmake -S . -B build
cmake --build build --target ducklips_contour_pressure_probe ducklips_contour_pressure_flow_probe
build/ducklips_contour_pressure_flow_probe
python tools/contour_pressure_flow_reference.py
ctest --test-dir build -R '^ducklips_contour_pressure_(storage|flow)$' --output-on-failure
```

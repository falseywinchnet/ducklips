# Continuous gas storage over the articulated contour

`ContourPressureStorage` implements a continuous gas field on the positive
air region of the existing medial contour. It replaces the *storage*
representation with two overlapping mass moments. Their sum is the total
air mass expressed as reference-density volume. It does not divide the air
at the old lower/upper chamber boundary.

A controlled compression sweep crosses that old boundary while retaining
both moments and total mass. The maximum relative density is below 1.314;
forces remain finite and integrated force work converges to the stored-energy
change. This removes the old partition singularity in this storage control.
It is not yet a new coupled source, a sustained-voicing result or speech audio.

## The same two shape functions describe area and gas effort

Let phi0,phi1 be the lower and upper contour shape functions: each is constant
in its end segment and linear between the two tissue centres, with
phi0+phi1=1. Signed area is a phi0+b phi1; the air domain contains only its
positive part. No gas is stored in occupied tissue.

The two gas coefficients u0,u1 describe a continuous enthalpy-effort field:

```
u(z) = u0 phi0(z) + u1 phi1(z)
r(z) = [1 + (gamma-1) u(z)/B]^(1/(gamma-1))
p(z) = B/gamma [r(z)^gamma - 1].
```

Here r is density relative to ambient, B is the ambient bulk modulus, and
p is pressure relative to ambient. The equation of state is unchanged.
The coefficients are not local pressures at two separated chamber centres.
A coefficient can be large outside the occupied air region without implying
that the air has that pressure; the physical field is evaluated on the air
region only. The bracket defining r must be positive on that entire region,
including its contact endpoint.

Define the pressure coenergy and the two mass moments by

```
Psi(u,a,b) = integral_air A(z) p(u(z)) dz
Mi = dPsi/dui = integral_air A(z) r(z) phii(z) dz
H(M,a,b) = u dot M - Psi.
```

H equals the integral of the existing isentropic available compression energy.
Because the shape functions sum to one, M0+M1 is total reference-mass volume.
M0 and M1 are overlapping moments: **they cannot be substituted for the old
chamber masses**. In particular, M0 remaining positive after contact crosses
the old boundary does not mean air remains below that boundary.

Differentiation of the same potential gives

```
dH = u dot dM + F0 da + F1 db
Fi = -integral_air p(z) phii(z) dz
Cij = dMi/duj = integral_air A r^(2-gamma)/B phii phij dz
Dij = dMi/daj = integral_air r phii phij dz.
```

The moving contact endpoint contributes no coenergy boundary term because
A=0 there. A whole zero-area plateau still has distinct one-sided shape
forces; the implementation chooses the occupied/closed side, consistent with
the contour convention. Away from that degeneracy, C is positive definite
when both shape functions remain independently represented in the air.

## Native implementation

The component evaluates moments, energy/coenergy, shape forces, compliance C
and the geometry derivative D. `recover` solves for the two effort coefficients
from prescribed moments using C and a bounded Newton line search. Failed
operations preserve the supplied output. Empty air support, invalid moments
and a nonpositive density bracket are rejected.

Integration splits at the actual contour knots and contact point. It uses
8/16-point Gauss comparison with bounded subdivision; it introduces no minimum
area or residual gas volume. Gas energy uses the existing cancellation-safe
`GasStorage` evaluation. The implementation uses fixed-size state and bounded
stack recursion; no heap storage is introduced by the component.

The API exposes storage and its derivatives, not a time advance. No production
renderer switches to it in this commit.

## Checks

[The independent report](contour_pressure_reference.json) evaluates 72 open,
lower-contact and upper-contact states, including both sides of the old
chamber boundary and much larger effort variations than the initial sweep.
It compares every returned field against separate 80-digit integration, then
recovers the coefficients from independently computed mass moments. The
maximum scaled field discrepancy is 1.84e-15, and recovered effort coefficients
differ by at most 4.09e-7 Pa. The report records exact source/binary/input hashes.

Native checks also verify the Legendre identity, positive compliance,
fixed-moment energy derivatives and output preservation on rejection.
The centred force-difference
check excludes a zero plateau, where a two-sided derivative is not defined.

In the mass-preserving sweep, upper signed area is 1e-6 m² and lower signed
area moves from -4.8e-6 to -5.2e-6 m². The old lower chamber disappears at
-5e-6 m². The new field begins at uniform 700 Pa enthalpy effort and retains
both of its initial moments throughout the sweep. This is a prescribed
geometry constitutive test, not a simulated tissue trajectory.

The report checks four subdivisions of the same path. Trapezoidal integration
of F0 da converges at second order to H_final-H_initial. The maximum density
is approximately 1.31303 times ambient, and there is no mass deletion or reset.
At 800 subdivisions, maximum moment error is 2.59e-25 m³ and the force-work
residual is 2.49e-15 J. Successive work-error ratios are essentially four.
These results validate this path and storage formulation; they do not establish
that an eventual coupled trajectory has physiologically realistic pressures.

[Validation receipts](contour_pressure_validation.json) record the focused
CTest, address/undefined-behavior sanitizer runs, and byte comparisons of
native versus sanitized fixture and sweep outputs.

## What the momentum representation must now respect

The mass variables changed meaning, so the old three-passage momentum system
cannot simply be attached to them. With fixed contour shape functions and
reference-mass flux Q(z), continuity gives

```
dM0/dt = Q_in - (1/ell) integral_between_centres Q(z) dz
dM1/dt = (1/ell) integral_between_centres Q(z) dz - Q_out,
```

where ell is the distance between tissue centres and Q is zero in occupied
regions. The internal transfer is therefore an integrated flux, not flux at
the old chamber seam. This follows by integration by parts of the local mass
continuity equation against each shape function.

The next coupled representation must derive its velocity/flux field, inertia
and losses consistently with that identity and with the same continuous gas
density. Moving-contact kinematics and complete loss of air support still
need treatment. The storage component supplies the required mass and force
relations, but does not claim to have solved those remaining dynamics.
Performance optimization remains parked; no pronunciation improvement is
attributed to these constitutive checks.

Reproduce with the repository's scientific Python environment:

```
cmake -S . -B build
cmake --build build --target ducklips_contour_pressure_probe
build/ducklips_contour_pressure_probe
python tools/contour_pressure_reference.py
ctest --test-dir build -R '^ducklips_contour_pressure_storage$' --output-on-failure
```

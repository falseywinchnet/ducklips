# Why the lower pocket stops the pressure-driven source

The source's 10.91 ms stop is not just a time-step problem. An independent
continuous solution of the same equations passes through a rebound, but the
lower pocket reaches roughly **4.29e11 times ambient density** and
**1.93e21 Pa** above ambient. These values reject this trajectory as a credible
speech-source result, even though its energy balance and numerical refinement
checks pass. This investigation does not add a speech recording or change
which NATO words are legible.

The likely defect is in the representation of the shrinking chamber and its
outflow: gas density is uniform within each chamber and jumps at their fixed
boundary. The viscous passage law uses those two densities all the way to that
boundary. Integrating a continuous pressure/density profile instead changes
how the residual air escapes. That isolates a consequential modeling assumption;
it does not yet validate a replacement dynamic source.

## Native continuation is insufficient

The new `ducklips_contour_pocket_probe` reproduces the 700 Pa trajectory and
saves the full tissue, gas, passage and network state around 10.8802 ms. On
failure it halves the time step; successful steps may grow back toward the
original step. This is a bounded diagnostic, with no error estimator. Its
`completed` status means it reached **11 ms**, not the earlier 100 ms target.

| Nominal steps per 100 ms | Diagnostic outcome | Minimum lower volume | Maximum relative density |
| --- | --- | --- | --- |
| 76,800 | Reached 11 ms after 14 rejected steps | 1.06e-23 m³ | 2.02e8 |
| 307,200 | Stopped when requested step fell below 1e-14 s | 3.63e-25 m³ | 2.78e9 |

The coarse run can be made to complete while producing a very different
compression peak. Completion and a small work residual therefore do not prove
trajectory accuracy. Neither diagnostic is promoted to a rendering algorithm.

## A local coordinate that resolves the pocket

Let a,b be the lower and upper signed openings. With the present thicknesses,
the opening at the fixed chamber boundary is S=(a+5b)/6. In the investigated
phase a<0, S>0 and b>0. The contact point is below the boundary. Let

```
k = (b-a)/0.0015
V_lower = S²/(2k)
z_lower_centroid = z_boundary - S/(3k).
```

The lower chamber disappears at S=0, although the upper chamber still contains
a substantial amount of air. This is the contact point crossing a numerical
partition, not both exterior openings closing.

The independent continuation uses log mass and log S, with exact local
integrals expressed through S. It avoids subtracting almost coincident axial
positions to resolve the wedge. In the checked region, the upper centroid is
in the constant upper-cap segment. For the middle passage's lower part:

```
I1_lower = log(3/2)/k
I3_lower = 5/(8 k S²).
```

The upper part contains `log(b/S)/k` in I1 and
`(1/S² - 1/b²)/(2k)` in I3, followed by the constant-cap contribution.
The script is restricted to this signed-contour phase and a symmetric body;
it is not a general replacement for the contour geometry implementation.

These local formulas and their mass/shape gradients agree with the separately
implemented clipped integrals at 80-digit precision in 120 cases, covering
S=1e-7 through 1e-16 m² and relative densities up to 1e12. Maximum scaled
error is below 9.33e-16. This validates the formula evaluation in the tested
regime, not the physical validity of such densities.

## Continuous equations also rebound

`contour_pocket_reference.py` integrates the existing coupled equations with
SciPy Radau from the saved checkpoint for 40 microseconds. It includes the
same gas potential, kinetic shape/mass derivatives, tissue/contact energy,
viscous and quadratic losses, and network coupling. The lower exterior cap
remains closed. The two symmetric tissue pairs are reduced to their area
coordinates without changing their mass, elasticity or damping.

The rebound is located by the zero crossing of dS/dt. At relative tolerances
1e-8 and 1e-9, the final scaled states differ by 3.28e-10; relative differences
in the rebound mass, volume, density and pressure are below 1.00e-6. An
independent run from the four-times-finer native checkpoint also reaches the
same extreme-compression regime. The maximum recorded continuous work
residual is 2.24e-15 J, against initial stored energy about 2.66e-5 J.

At the fine rebound, lower volume is about 4.33e-28 m³ while its retained
reference-density mass volume is still 1.86e-16 m³. The residual air creates
a large elastic gas-energy barrier. That mechanism exists in the continuous
model; it is not removed by choosing a more accurate native time integrator.
These are mathematical model outputs, not predicted physiological pressures.

## Testing the density assumption in the viscous loss

The existing loss coefficient is

```
R_old = 12 eta w² (I3_lower/r_lower² + I3_upper/r_upper²),
```

where eta is dynamic viscosity, w is fold width, and r is density relative
to ambient. Near disappearance, both I3 terms diverge as S^-2. The upper
chamber remains near ambient density even while the lower one becomes very
dense. Treating that whole upper half-passage as ambient gas can impose a
large resistance immediately next to the narrow seam.

A useful steady constitutive reference follows from the parallel-plate
Poiseuille relation and the existing isentropic equation of state. The
parallel-plate starting point is covered by the [MIT Couette/Poiseuille
notes](https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_couet_and_pois/).
The compressible integration below is derived here. For reference-mass flow
Q=r q, where q is actual volume flow:

```
dp/dz = -12 eta w² Q/(r A³)
p(r) = B/gamma (r^gamma - 1)
mu(r) = B/(gamma-1) (r^(gamma-1) - 1)
Phi(r) = B/(gamma+1) r^(gamma+1)

Q = (Phi(r0)-Phi(r1)) / (12 eta w² I3)
R_pressure = 12 eta w² I3 (mu(r0)-mu(r1))/(Phi(r0)-Phi(r1)).
```

Here mu is the gas mass effort and B the ambient bulk modulus. Phi is obtained
by integrating r dp. The coefficient is positive, symmetric under swapping
endpoints, and reduces to `12 eta w² I3/r²` at equal density. Its loss power
is nonnegative. The implementation evaluates the ratio with logarithms and
`expm1` to retain the equal-density limit. 125 cases agree with 80-digit
power differences to relative error below 4.38e-16.

A diagnostic replaces only the viscous coefficient in the continuous solve
with this pressure-profile value. Tissue parameters, kinetic geometry and
other forces remain the same. It stops at specified positive seam openings;
those are observation endpoints, not an area floor or a state reset.

| Observed seam area | Lower reference-mass volume | Lower relative density | Upper area rate |
| --- | --- | --- | --- |
| About 1e-13 m² | 3.64e-21 m³ | 1.06e3 | -0.00710 m²/s |
| About 1e-15 m² | 9.27e-24 m³ | 2.70e4 | -0.00710 m²/s |

The residual mass keeps decreasing and the earlier rebound is absent up to
these observation endpoints. However, density and pressure still grow to
extreme values. This does not establish a valid zero-volume limit, complete
closure, or realistic phonation. The pressure profile is an exact *steady*
viscous relation under the declared slit/EOS assumptions; using it during
unsteady motion is itself an approximation. Passage inertia still uses the
old piecewise-density representation. The exterior tract remains the existing
linear acoustic load.

## Consequence for the next engine change

The next work is the spatial gas representation and disappearing-cell
transport, not further tuning of contact damping or reducing the time step
until this particular bounce renders. In particular:

- Pressure/density variation inside the narrow passage must be represented
  consistently in both inertia and dissipation.
- A numerical chamber boundary must be allowed to move or its state must be
  transported when contact sweeps it away. Its mass and energy cannot simply
  be deleted.
- A complete trajectory must converge in motion and compression, in addition
  to conserving mass/work, before using it for a new speech sample.

This is a diagnosis and a checked constitutive alternative, not an implemented
replacement source. The production source and latest spatial-pressure audio
are unchanged. Performance optimization remains parked.

[The numerical report](contour_pocket_reference.json) contains the independent
checks and continuous results. [Validation receipts](contour_pocket_validation.json)
record source, binary and output hashes, native failures and sanitizer checks.
The sanitized native diagnostic is byte-identical to the release result.

Reproduction, using the repository's scientific Python environment:

```
cmake -S . -B build
cmake --build build --target ducklips_contour_pocket_probe
mkdir -p out/contour-pocket
build/ducklips_contour_pocket_probe 76800 > out/contour-pocket/native-76800.jsonl
build/ducklips_contour_pocket_probe 307200 > out/contour-pocket/native-307200.jsonl
python tools/contour_pocket_reference.py --out out/contour-pocket/continuous.json
python tools/contour_pocket_reference.py --rtol 1e-9 --out out/contour-pocket/continuous-fine.json
python tools/contour_pocket_reference.py --input out/contour-pocket/native-307200.jsonl --rtol 1e-9 --out out/contour-pocket/continuous-finer-start.json
python tools/contour_pocket_reference.py --pressure-loss --rtol 1e-9 --out out/contour-pocket/pressure-loss.json
python tools/contour_pocket_reference.py --pressure-loss --stop-seam 1e-15 --rtol 1e-9 --out out/contour-pocket/pressure-loss-small.json
python tools/contour_pocket_audit.py
```

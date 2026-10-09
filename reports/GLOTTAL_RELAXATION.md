# Resolve viscous decay without damping every tissue mode

The continuous source now has an explicit `relaxation` update which fits the
linear viscous decay of its coupled fluid/tissue momentum modes. It retains
the shared inertia and physical force law. It is not the default update or a
new speech renderer. The experiments below verify this component and also show
that it does **not** by itself resolve the coupled contact-trajectory error.

## Separate contact crossing from the approach to contact

New controls start from the same complete physical source and graph state,
saved at 8 ms of the previous 19,200-step tube run, and advance another 16 ms.
This removes differences accumulated before that point from the comparison.
`continuous_tract_probe --graph-state` accepts nine symmetric source variables
followed by eight graph charges and eight branch momenta.

The optional `--contact-step` limits the dissipative crossing step. A failed
fixed-topology step is subdivided until a crossing can be attempted below that
limit; every substep solves a new graph response. It does not floor an area,
erase stored gas, or alter a force. The experiments use a maximum crossing
step of approximately 78.125 ns, preserving the requested maximum step for
surrounding motion.

This reduces crossing loss but does not establish accuracy. Localized controls
with maximum open steps of 2.5 and 1.25 microseconds still differ by 210.10 Pa
in final tract pressure and by 5.55 in the declared source-state normalization.
The unlocalized 2.5 microsecond control differs from the localized fine control
by 109.48 Pa. Localizing the crossing alone is not a successful trajectory repair.

## The open pressure equations become extremely stiff

Independent local Jacobians use the continuous-field equations at four saved
open source states, holding tract effort at its recorded value. This is a local
source-timescale diagnostic, not equilibrium stability or a full coupled-graph
eigensystem. The [contact receipt](continuous_tract_contact_reference.json)
records quadrature and difference-step refinements.

| Time in original run | Lower area | Most negative local real part |
|---|---:|---:|
| 8.000 ms | 1.44726 mm2 | -4.45618e4 /s, complex mode |
| 8.300 ms | 0.174230 mm2 | -1.05983e6 /s, real mode |
| 8.400 ms | 0.0157394 mm2 | -3.25208e7 /s, real mode |
| 8.425 ms | 0.00169948 mm2 | -4.51964e8 /s, real mode |

The last decay time is about 2.21 ns. For scalar decay at that rate, the
symmetric midpoint factor at a 1.25 microsecond step is about -0.992945,
whereas continuous decay is practically zero. The symmetric method can
preserve an alternating numerical remnant of a mode that should decay within
the step. Subdividing only the contact crossing cannot fix a stiff mode already
present in the open phase. Other local modes still oscillate near 100 kHz;
fitting viscous decay does not make those phases exact.

## Fit resistance modes in the shared inertia metric

Let `A` be the average endpoint inverse inertia and `R` the midpoint linear
resistance in normalized coordinates. Factor `A = L L^T` on active coordinates
and diagonalize `S = L^T R L`. Apply the checked scalar function to each rate:

```
beta(x) = coth(x/2) - 2/x,     x = step * rate
B = L beta(step * S) L^T
v_average = (A + B) p1 / 2 + (A - B) p0 / 2
```

For constant inertia and linear resistance this reproduces exact momentum
decay, including integrated velocity. Low rates have `beta(x) = x/6 + O(x^3)`;
stiff rates tend to one. The weighting acts on coupled modes, rather than
imposing a global weight from the fastest mode on all tissue motion.

The implementation uses active-coordinate Cholesky factorization and a
seven-dimensional symmetric Jacobi eigensolve. Eliminated coordinates remain
exactly zero. Symmetry checks use the diagonal quadratic-form scale: an
off-diagonal can cross zero by cancellation of signed products, so its nearly
zero result is not its floating-point uncertainty scale. Material asymmetry
and negative resistance are rejected; negative eigenvalues within the stated
arithmetic uncertainty are treated as zero.

The symmetric finite configuration-energy derivative is retained. With
`Delta p = p1 - p0`, the additional nonnegative step-work term is

```
D_bias = Delta p^T B Delta p / 2.
```

This preserves the source work identity up to the measured solver/arithmetic
residual. The term uses the existing `numerical_loss` field. In constant-inertia
linear decay, where the step is exact, it equals viscous dissipation missed
by replacing varying velocity with its step average. In a nonlinear moving
source it is a reported work defect, not proof that unresolved physical heat
has been calculated exactly.

The fitted update preserves active port constraints. Changes of constraint
still use the existing dissipative transition and its loss accounting. Graph
port work continues to cancel to floating-point accuracy. No waveform is
prescribed.

## Passing checks and the accuracy result that did not pass

The [independent receipt](glottal_relaxation_reference.json) contains 80 seeded
coupled inertia/resistance cases, including eliminated coordinates. Native
matrix functions agree with independent eigendecomposition within 1.36e-13
in the inertia normalization. The resulting linear transfer agrees with an
independent matrix exponential within 1.11e-12. The tested bias matrices lie
between zero and inverse inertia.

Fresh nonlinear coupled runs use 10, 5 and 2.5 microsecond steps over 2 ms,
with fixed and moving tract metrics. Endpoint refinement ratios are
3.972/4.027 and 4.282/4.072. They agree with independent continuous equations
under refinement. The ODE artifacts are hash-verified results from the previous
audit; the native fitted runs are new. Coarse whole-trajectory startup errors
remain substantial: viscous fitting does not resolve ultrasonic pressure
oscillations by itself.

The common-state contact experiment gives:

| Fitted maximum step | Port closures / reopenings | Bias plus crossing work defect |
|---|---:|---:|
| 5 microseconds | 3 / 3 | 9.14654e-8 J |
| 2.5 microseconds | 3 / 3 | 1.54323e-8 J |
| 1.25 microseconds | 3 / 3 | 4.14947e-9 J |

All runs complete and satisfy work/mass checks. Nevertheless, adjacent fitted
refinements differ by 90.87 and 206.70 Pa in final tract pressure, and by 8.65
and 5.06 in normalized source state. This is not a converged sequence. On the
finest grid, localized symmetric and fitted updates differ by 3.83 Pa and
0.0723 normalized source units. Agreement between methods at one grid does
not establish accuracy when their own refinements remain large.

Ten focused CTests pass in Release and AddressSanitizer/UndefinedBehaviorSanitizer.
New checks cover fitted coupled work, sealed ports, allocation-free steps,
rotated modes, eliminated coordinates and cancellation-scale symmetry.
Existing source, contact, onset and throughflow controls pass as well.

The remaining problem concerns full fast pressure evolution and its interaction
with narrowing geometry and the tract. This update addresses one verified
decay limit; it is not a reason yet to replace the speech path or claim better
NATO intelligibility. The tube also remains inviscid with an ideal pressure
outlet, so it is not a complete human tract. Performance optimization is parked.

## Reproduce

```sh
cmake -S . -B build
cmake --build build --target ducklips_glottal_relaxation_probe ducklips_graph_continuous_glottis_probe ducklips_continuous_tract_probe
ctest --test-dir build --output-on-failure -R 'glottal_relaxation|graph_continuous_relaxation'
OPENBLAS_NUM_THREADS=1 python tools/glottal_relaxation_reference.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_tract_contact_reference.py
```

The contact receipt stores its physical checkpoint and diagnostic states, so
a clean checkout can reconstruct inputs without old local trajectory files.
An archived control executable is used for the pre-fitting controls when
present; otherwise the current symmetric path is used. Executed commands,
binary/source hashes and output hashes are recorded.

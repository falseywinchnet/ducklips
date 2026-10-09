# Constrained tissue experiment: not yet an engine replacement

The rigid-body issue in [TONGUE_GEOMETRY.md](TONGUE_GEOMETRY.md) cannot be
resolved by declaring a new contour volume-preserving. This experiment tests
that premise against the native Quebec starting pose. It finds that a model
can preserve the **sum** of tissue areas while substantially inflating local
pieces. The stricter elementwise-conservation solve fails near closure.
Neither candidate is installed in the speech renderer.

## Physical question and explicit assumptions

A replacement should form a dorsal closure through tissue deformation,
respect the palate, and account for tissue redistribution. The biological
motivation is the tongue's muscular-hydrostat behavior. Hyperelastic tissue,
near incompressibility, muscle forces and contact with surrounding structures
are established ingredients of biomechanical tongue models; see Calka et al.,
[Biomechanical modeling of the human tongue](https://arxiv.org/pdf/2407.13259),
sections 4.3–4.4. That source motivates the constraints; this experiment does
not reproduce its anatomy, constitutive calibration or muscle model.

The native probe exports the outer wall, target, rest pose and current rigid
poses. The measured centre is (3.059093051, 5.322166367) cm at 1.310 s in the
previous Papa/Quebec diagnostic render. The target lies 6 cm from the lips,
with zero gap. The experimental patch spans 0.9 cm horizontally, using the
existing K gesture's width value as a fixture assumption. This interpretation
is different from the old Gaussian-width parameter. Patch endpoints are
identical across mesh refinements.

The candidate uses a triangular sagittal slice of the current circle. Lower
nodes, below centre height minus half the radius, stay fixed; the dorsal
patch is prescribed toward the palate. Those anchors and prescribed surface
motion are **not measured anatomical insertions or muscle actuation**.
Out-of-plane thickness is constant. The blade, tip, lips, fitted vowel
residuals and acoustic state are absent from this isolated experiment.

For each triangle, `F` is the deformation gradient and `J = det(F)`. The
normalized plane-strain energy is

```
W = (tr(F^T F) - 2)/2 - log(J) + lambda*log(J)^2/2
lambda = 2*nu/(1 - 2*nu), shear modulus = 1
```

The absolute shear modulus cancels in this displacement-controlled static
problem. No absolute force or dynamic timing is inferred. The default
Poisson ratio 0.49 is a near-incompressible modeling choice, not a fit to the
NATO speaker. The comparison at 0.499 tests sensitivity, not a new calibration.

Two equality-constraint choices are implemented: preserve the summed slice
area, or preserve each triangle's area. Both impose unilateral outer-wall
contact. Inverted trials, violated constraints and unconverged solves are
rejected and recorded, never emitted as a usable pose. The last valid mesh
is saved separately when continuation fails.

## Results and rejection

The globally constrained model reaches the prescribed closure with zero
reported wall penetration and fixed anchors on all three meshes. Across the
33 solves, the worst relative total-area error is 9.42e-11. However:

| Nodes | Minimum local J at closure | Maximum local J at closure | Isolated body airway loss, mL | Rigid comparison loss, mL |
|---:|---:|---:|---:|---:|
| 51 | 0.9246 | 1.1351 | 4.364 | 29.394 |
| 123 | 0.9454 | 1.6805 | 4.269 | 29.329 |
| 227 | 0.9205 | 1.4117 | 4.231 | 29.307 |

`J = 1.6805` means **68.05% local area expansion**, compensated elsewhere.
Thus global volume preservation is not a sufficient tissue-conservation
check. The apparent airway improvement is not accepted as a physical speech
improvement. With nu = 0.499 on the 123-node mesh, maximum J falls to 1.2140,
but it still represents 21.40% local expansion. Simply increasing the bulk
penalty does not establish the desired incompressible limit.

With each triangle constrained individually, the 51-node model is valid
through activation 0.8 using steps of 0.1; the next solve fails. Refining
continuation to 0.01 reaches activation 0.88, with local J within 1.7e-13 of
one and dorsal gaps approximately 0.78–0.85 mm. The solve at 0.89 fails and
contains an inverted trial element. The 123-node model also fails at 0.9 with
0.1 steps. These are **solver/discretization failures**, not evidence that a
human tongue cannot make the closure. Constant-strain triangles with one
incompressibility constraint each leave few deformation degrees of freedom;
locking is a concern, but its responsibility is not proved by these runs.

The native rigid circle additionally penetrates the sampled outer wall by
about 0.49 mm at closure in this fixture. The existing airway code turns
contact into zero gap; it does not redistribute tissue to satisfy the wall.

Airway losses in the table use the same isolated-body ray intersection and
native gap-to-area power law for both candidates. They are **not** the full
oral-volume numbers from the speech trace. More fundamentally, conservation
of sagittal tissue area does not establish conservation of 3D tissue volume
under that empirical nonlinear gap-to-area law. A production replacement
needs an explicit transverse geometry or a justified volume relationship.

## Verification and reproducibility

Analytic energy, global-area, element-area and contact derivatives are checked
against directional finite differences. The largest absolute discrepancy in
the three-mesh run is 4.04e-9. Projected stationarity residuals are checked after
successful solves. Ray tests cover an outside hit, a miss and a wall origin
inside the pressed body; an exit intersection must not be reported as air.
The Python wall is checked against the native fixture. The native exporter
was compiled with ASan/UBSan and its output matched the release fixture byte
for byte. No speech audio or spectrogram is regenerated for this experiment.

```sh
cmake --build build --target ducklips_tongue_geometry_probe
./build/ducklips_tongue_geometry_probe 3.059093051 5.322166367 6 0 \
  > out/tongue_geometry/tissue-fixture.json
OPENBLAS_NUM_THREADS=1 python tools/tongue_tissue_reference.py \
  out/tongue_geometry/tissue-fixture.json --out out/tongue_tissue --mesh 7 11 15
OPENBLAS_NUM_THREADS=1 python tools/tongue_tissue_reference.py \
  out/tongue_geometry/tissue-fixture.json --out out/tongue_tissue_stiff \
  --mesh 11 --poisson .499
# These runs intentionally return 1 when an invalid continuation trial is rejected:
OPENBLAS_NUM_THREADS=1 python tools/tongue_tissue_reference.py \
  out/tongue_geometry/tissue-fixture.json --out out/tongue_tissue_local \
  --mesh 7 --incompressible
OPENBLAS_NUM_THREADS=1 python tools/tongue_tissue_reference.py \
  out/tongue_geometry/tissue-fixture.json --out out/tongue_tissue_local_fine \
  --mesh 7 --incompressible --steps 100
```

The [receipt](tongue_tissue_reference.json) embeds results and hashes for all
five cases, including the 123-node local-constraint failure. It preserves the
native fixture as well. Geometry figures and meshes remain regenerable under
`out/`. Performance work is still parked. Next geometry work must resolve local
conservation and contact through closure, then connect that tissue geometry to
three-dimensional airway volume before producing speech-level claims.

The [continuous native material-map follow-up](TONGUE_FLOW.md) now reaches
this contact patch with local conservation and the same lower anchor. It
establishes a conservative kinematic solution despite the triangle-solver
failure; calibrated tissue mechanics and speech integration remain open.

# Native conservative material map

The previous [triangle experiment](TONGUE_TISSUE.md) failed near closure when
it constrained every element's area. A continuous material map now reaches
the same Quebec contact patch while keeping the lower tissue fixed and
preserving local volume. Therefore that failure does not establish that the
fixed-width closure itself is impossible. It remains a failure of that
particular discretization/solve; these results do not isolate locking from
other solver limitations.

`ConservativeTongueMap` is a native C++ component. This report records its
isolated validation. The subsequent [speech integration](CONSERVATIVE_TONGUE_SPEECH.md)
connects it to `TractModel`, including the blade attachment and airway
intersections. It supplies conservative kinematics, not calibrated tissue
forces or a human pronunciation result.

## How conservation is obtained

Let `f(x)` be a commanded dorsal lift, `F(x)` its antiderivative, and `S(y)` a
smooth gate that is zero in anchored tissue and one above the moving region.
The material velocity, with respect to a dimensionless command interval, is

```
vx = -F(x) S'(y)
vy =  f(x) S(y)
```

Its divergence is identically zero. Raising the dorsum is coupled to
redistribution of the surrounding tissue. There is no after-the-fact scale
factor to repair lost or created area. Extending the map as `(X,Y,z)` also
preserves three-dimensional volume locally, although it supplies no
transverse motion or measured transverse anatomy.

The profile uses cubic Hermite pieces, with their exact polynomial
antiderivatives. The gate is a quintic smoothstep. Implicit midpoint advances
material points. For its trace-zero local derivative `A`, the derivative of
a step is

```
M = (I - h A/2)^-1 (I + h A/2)
```

In two dimensions the numerator and denominator have equal determinants, so
`det(M) = 1`. The Newton solve is checked for convergence and finite values.
The native map also propagates the material deformation gradient, which is
checked independently by perturbing input positions. Reversing the command
with the same step count reverses the discrete material map.

The map uses fixed storage and performs no dynamic allocation. The command
parameter is not physiological time. No runtime throughput claim is made.
The caller must supply a contact-compatible profile; the general mapping API
does not itself enforce arbitrary palate, jaw or neighboring-tissue contacts.

## Quebec fixture

The initial body and wall come from the native exporter used in the previous
experiment. A 0.9 cm dorsal patch is driven to the palate, with smooth 1.5 cm
shoulders. The lower anchor remains at `centre_y - radius/2`, as before. The
gate reaches one at the lowest resting height of that patch. This avoids
forcing the deformation into an unnecessarily thin transition layer.
These profile and anchor choices are explicit geometric assumptions, not
anatomical or muscle measurements.

At full command, on 12,417 material and contact samples:

| Check | Result |
|---|---:|
| Maximum local volume error, `abs(det(F)-1)` | 2.35e-14 |
| Independent finite-difference volume error | 4.43e-9 |
| Independent gradient absolute error | 4.00e-8 |
| Forward/inverse position error | 2.58e-14 cm |
| Anchor displacement | exactly zero |
| Maximum sampled wall penetration | 3.02e-14 cm |
| Maximum absolute patch clearance | 3.30e-11 cm |
| Maximum principal stretch | 3.045 |

The test covers eleven command levels through full contact. Halving the
integration step reduces the largest position difference by approximately
four, consistent with second-order integration: 16→32 steps gives
7.45e-4 cm, 32→64 gives 1.86e-4 cm, and 64→128 gives 4.66e-5 cm.

Drawing the mapped surface with straight polygon edges introduces a separate
geometric approximation. Its relative area change decreases from 1.03e-5
with 512 vertices to 6.56e-7 with 2,048 and 4.11e-8 with 8,192. This must not
be confused with the pointwise conservation of the material map. A later
engine integration must control both errors rather than infer conservation
from a plot or coarse polygon alone.

The isolated-body ray/gap calculation gives 72.196 mL at rest, 63.369 mL for
this mapped body, and 42.911 mL for the rigid translation. Those are diagnostic
outputs of the existing empirical gap-to-area law, not validated physical air
volumes or evidence of improved sound. The three-dimensional airway relation
identified in the previous report still needs to be resolved.

## Limits and integration requirements

A principal stretch of about three is not justified by the conservation
identity. This map has no calibrated elastic energy, muscle stress, tissue
inertia or physiological time law. Its contact and conservation success is a
kinematic existence result. It does not establish that this is the motion a
human tongue would choose or could sustain under a particular command.

Integration requires consistent attachment of the blade and tip to the
deformed body, geometric treatment of combined gestures, and continuous
intersection of the mapped surface with airway sections. Connecting only the
new upper contour while leaving the blade attached to the old translated
circle would produce a disconnected articulator and is not acceptable.
This isolated component did not supply a speech sample. The subsequent
[integration report](CONSERVATIVE_TONGUE_SPEECH.md) records actual speech from
the connected engine and its remaining reach limitation. Performance
optimization remains deferred by user instruction.

## Reproduce and validate

```sh
cmake --build build --target ducklips_tongue_flow_probe ducklips_tongue_flow_tests
ctest --test-dir build -R tongue_flow --output-on-failure
OPENBLAS_NUM_THREADS=1 python tools/tongue_flow_reference.py \
  out/tongue_geometry/tissue-fixture.json --out out/tongue_flow
```

`tools/tongue_flow_reference.py` generates the contact profile and checks the
native executable's returned coordinates. It saves full material samples,
deformation gradients, a full geometry figure, hashes and a numerical
receipt. The native regression also checks fixed tissue, exact upper-region
shear, reverse mapping and invalid configuration/command handling.

The Release suite passed all eight tests. The final native flow regression
also passed under ASan/UBSan. The exact measurements and source hashes are in
[tongue_flow_reference.json](tongue_flow_reference.json).

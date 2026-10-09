# Near-closure compression and retained glottal air

The rectangular chamber model's strong compression persists under time-step
refinement. A new `GlottalContour` geometry reference separates local contact
from air storage: an outlet can seal exactly while a wedge of air remains
beside the sloping medial surface. It supplies integrated storage, pressure
force areas and passage metrics from that one surface.

The pressure and motion measurements below come from the rectangular
chamber engine. The subsequent [open-contour dynamics](CONTOUR_GLOTTAL_DYNAMICS.md)
couple the new geometry in `GraphContourGlottal`; contact transport remains
unfinished. Neither reference has an audio or pronunciation result attached
to it.

## What is causing the compression in the current equations?

The motion probe now captures both chambers' signed openings and each
chamber's peak compression at every integration step. It also records actual
net mass transfer and volume change, rather than relying on its sparse CSV
samples to catch short pressure peaks.

With unchanged 700 Pa reservoir pressure and body parameters:

| Steps over 100 ms | Maximum density / ambient density |
| --- | ---: |
| 307,200 | 2.25129 |
| 1,000,000 | 2.25310 |

The upper chamber's million-step peak occurs at 61.35135 ms. Its volume is
`1.11653e-11 m³`, gas pressure is approximately **211 kPa relative to ambient**,
and its volume is shrinking at `6.098e-7 m³/s`. Air is leaving both ends:
`6.875e-7 m³/s` towards the lower chamber and `6.966e-7 m³/s` towards the
tract, measured in reference-density volume units. This is not an exactly
sealed cavity or a stale nonzero momentum behind a closed port.

The current slit resistance grows with inverse aperture cubed, including
the density conversion. At that peak, differential resistances on the two
escape paths are about `1.95e11` and `1.94e11 Pa s/m³`. A local gas-only
linear estimate with fixed neighboring efforts gives a pressure-relaxation
time of **12.7 µs**, while `V/|dV/dt|` is **18.3 µs**. This estimate omits
passage inertia, kinetic-head derivatives and neighboring motion; it is
not a measured eigenmode. It nevertheless shows that escape is not fast
relative to the chamber's contraction in these equations.

The lower chamber peaks earlier, at 60.45435 ms, with density ratio 1.74415
and relative gas pressure about 118 kPa. Its analogous relaxation/contraction
times are 21.0/40.8 µs. These are gas pressures in a model without tissue
collision. They must not be presented as tissue-impact measurements.
Experimental intraglottal recordings distinguish impact, preopening and
open-flow phases; see [Jiang and Titze's pressure measurements](https://pubmed.ncbi.nlm.nih.gov/8061769/).

The evidence establishes that the large compression is a feature of the
current volume/flow equations, rather than disappearing with finer stepping.
It does not establish that either the rectangular chamber geometry or its
quasisteady loss law is physiologically correct. No extra damping or aperture
floor has been added to conceal it.

## One explicit medial contour

The reference uses the existing lower/upper tissue-centre openings `a,b` and
longitudinal thicknesses `h0,h1`. The signed area is constant from each cap
to its tissue centre and affine between the centres:

```
z:       0          h0/2             h0+h1/2          h0+h1
area:    a -------- a --------------- b -------------- b
actual air area: max(signed area, 0)
```

This is a declared profile assumption based on the two tissue coordinates,
not measured anatomy and not the old soft minimum relabeled as geometry.
The fixed chamber boundary is at `z=h0`; air centroids are computed from the
positive-area integrals and move with the contour. The geometry exposes:

- Actual volume and air centroid in each chamber, including exact emptiness.
- Derivatives of each volume with respect to both signed control areas.
  Multiplying by glottal width gives the corresponding tissue force areas.
- Inlet, chamber-seam and outlet apertures.
- Integrals of `1/A` and `1/A³` along cap-to-centroid and centroid-to-centroid
  paths, separately attributed to each chamber's density.

Contact is the nonpositive part of the surface. A path reaching an exact
zero opening is sealed, even when adjoining storage is positive. No minimum
area or volume is introduced. If a cell becomes completely occupied, its
volume and adjoining passage metrics become zero; transporting its previous
physical contents remains the dynamic owner's responsibility.

For example, `h0=2.5 mm`, `h1=0.5 mm`, `a=3 mm²`, `b=0` seals the outlet but
retains **0.0625 µL** of air in the upper wedge, connected to the lower
chamber. The current rectangular rule would give that upper chamber zero
volume at the same local upper opening.

At the exact recorded peak-compression poses, the comparison is:

| Peak's chamber | Rectangular volume | Contour volume at the same two openings |
| --- | ---: | ---: |
| Lower | 0.2588 µL | 5.2961 µL |
| Upper | 0.011165 µL | 0.017345 µL |

These are alternative geometric measurements of the same pose. They are
not a corrected pressure trajectory. Retrofitting these volumes while
keeping old air masses or forces would change energy without accounting
for the change.

## Verification and numerical contact handling

The [independent quadrature receipt](glottal_contour_reference.json) covers
**1,131 valid profiles**, including 1,001 poses from the million-step motion
control, signed random profiles and explicit contact/empty-cell cases.
Two invalid geometries reject. Scipy quadrature operates on an independently
interpolated area field, split at contact roots; it does not use the native
clipped-trapezoid formulas. Thin reciprocal-area boundary layers receive
geometric sampling subdivisions before numerical integration.

Maximum normalized volume error is `3.34e-16`, centroid error `4.45e-16`,
and force-area derivative error `1.26e-16`. Relative errors in the `1/A` and
`1/A³` path integrals are below `5.54e-10` across the tested profiles.
The scales and exact source/input/output hashes are in the receipt.

The gradient audit caught a real contact-edge error: computing a plateau's
endpoint weight through a rounded ratio could invent a tiny positive air
area there and assign a finite pressure-force area to the whole closed
plateau. The kernel now handles the constant end regions explicitly. At an
exact zero-area plateau, the volume derivative is the closed-side value;
a unique two-sided derivative does not exist at that event. This convention
is documented, not hidden behind a smoothing width.

Native contact/empty-cell tests and the complete profile batch pass
AddressSanitizer/UndefinedBehaviorSanitizer. Eight focused CTests pass.
[The compression/validation receipt](glottal_contact_validation.json) records
current diagnostic runs and the peak-pose comparison provenance.

## Next coupling requirement

The shape changes more than volume. In particular, the lower chamber's
pressure can act through the upper tissue coordinate and vice versa.
For density ratio `r_i=M_i/V_i`, the open-path metrics imply

```
L_j = rho0 sum_i [integral_ji(dz/A) / r_i]
R_j = 12 viscosity width² sum_i [integral_ji(dz/A³) / r_i²]
K_j = lambda_j² / (2 L_j)
```

The moving centroids and shape integrals make kinetic energy depend
explicitly on tissue coordinates. Its derivatives and the gas volume
forces must therefore enter the same work balance before this geometry
can replace the rectangular model. Closed-path momentum handling and mass
transport at full cell disappearance also remain required. The geometry
alone does not resolve jet separation, collision or physiological realism.

## Reproduction

```sh
cmake --build build --target ducklips_glottal_geometry_motion_probe ducklips_glottal_contour_probe
mkdir -p out/glottal-contact
build/ducklips_glottal_geometry_motion_probe 307200 out/glottal-contact/motion-307200.csv
build/ducklips_glottal_geometry_motion_probe 1000000 out/glottal-contact/motion-1000000.csv
python tools/glottal_contour_reference.py
```

The motion tool prints its full-step peak diagnostics as JSON on stdout and
writes sparse overview samples to CSV. The CSV must not be used to estimate
the narrow pressure maxima. No ordinary or enhanced spectrogram is generated
because no new speech audio was rendered. Performance optimization stays
parked.

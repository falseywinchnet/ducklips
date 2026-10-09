# Connect the reconstructed mouth to the articulated lip tube

The native moving-airspace reference now connects the front-mouth regions to
one shared lip tube whose section areas and lengths come from the existing
articulator replay. Its narrowing displaces air, its volumes store pressure,
and its air columns retain momentum. Explicit replacement seams remove the
lip and epilarynx lengths from the planar core before attaching a tube.

This is an inviscid integration reference over three separately initialized
1 ms Foxtrot intervals. It has not rendered new speech. The latest listening
evidence remains November improved and Foxtrot jammed around /tr/ in the
spatial-pressure recording. Complete lip closure, laryngeal coupling, losses
and continuous full-utterance state transport remain unfinished.

## Use the body model's partition and final section geometry

`ducklips_airspace_geometry_probe --core BODY` obtains the existing section
partition from `TractModel`, trims the domain at the corresponding reference
arc positions, and records those boundaries in the pose JSON. The reference
body has 64 oral sections, eight epilaryngeal sections and four lip sections.
Its planar core runs from 2.1875 cm to 16.40625 cm along the reference midline.
The core plus those two end segments accounts for the resting 17.5 cm exactly
once. Protrusion adds length through the lip sections.

The exporter scales coordinates, grid and extrusion depth by the body's
existing length scale. A 1.1 scale control verifies linear coordinate,
quadratic aperture and cubic volume scaling. This is the model's uniform
body-size scaling, not a fit to an individual person's anatomy. The original
22 mm reference extrusion and the diagnostic tongue/blade construction remain
model assumptions.

At all 66 saved poses, independent polygon clipping gives exactly the same
trimmed domain, with unchanged solid contours. The untrimmed export mode
reproduces the historical poses exactly. The coarse acoustic grid now ends at
the mouth seam, so its final cell centre stays behind the connection.

The lip data are the **final assembled** areas and lengths in the existing
native replay, including the lip-slot constriction. They are not inferred
from the mouth cap or replaced by a single nominal `lip_area`. For example:

| Foxtrot time | Lip section areas, from mouth to exterior |
|---|---|
| 5.084 s | 587.312, 522.449, 373.500, 283.611 mm² |
| 5.085 s | 586.554, 517.818, 359.975, 264.718 mm² |

Those four sections lose 0.103379 mL of volume over this interval. Each
section's volume change contributes to its wall-displacement flow. Section
areas and lengths are interpolated separately between the saved endpoints;
volume is their product. This is a declared interpolation of the native
control trajectory, not a claim that the articulator dynamics were evaluated
at every acoustic stage.

## Shared pressure and air-column inertia

`ducklips_planar_inertive_motion_probe --lip-tube` adds a pressure/storage node
for each lip section and uses `MovingGraphReference`. Every front-mouth cap
attachment connects to the first shared lip node. Internal mouth ports retain
their individual region pairs even where more than one pair uses a grid face.

For a lip section with area `A` and length `L`:

```
V = A L
C = V / (rho c²)
inverse_inertance_between_sections = 1 / [rho (L_left/2A_left + L_right/2A_right)]
wall_flow = -(V_end - V_start) / dt
```

The mouth-to-lip branch includes both the core path to the seam and the first
half of the lip section. That half-section is divided among incoming branches
in proportion to their exposed cap areas. Giving every branch the full lip
area would duplicate the tube's flow capacity. A symmetric split/unsplit
control verifies that dividing the core and cap across width leaves the shared
pressure response unchanged, including with moving lip area and length.
This allocation is a coarse junction model; it is not a resolved three-
dimensional mouth-to-lip flow field.

The terminal half-section connects to ideal ambient pressure. There is no
radiation resistance, end correction, turbulence or viscothermal loss here.
The epilarynx has been removed from the planar core but is not attached in
these controls: the new inlet seam is sealed. Each case starts with the
declared rear-core pressures of 600 Pa and front/lip pressures of 0 Pa; the
three cases are not a continuous evolution from the utterance's beginning.

Positive lip areas and lengths are required throughout the interval.
Complete closure is rejected before stepping. No positive aperture or volume
floor substitutes for missing zero-storage transport.

## Geometry, conservation and temporal resolution

The [receipt](planar_lip_reference.json) records independent audits of closure
(5.031–5.032 s), release (5.062–5.063 s), and lip narrowing (5.084–5.085 s).
Twelve runs at 16, 32, 64 and 128 subdivisions contain 728 steps and 2,196
endpoint/stage geometries. Independent Boolean geometry checks volumes,
internal openings and cap ownership. Geometric continuation witnesses retain
all five core region identities. A separate full charge-and-momentum
collocation system carries its own state from the initial conditions using
the independent metrics.

| Maximum discrepancy across those runs | Value |
|---|---:|
| Pressure against the independent evolution | 3.21e-9 Pa |
| Branch momentum | 1.19e-13 Pa s |
| Average branch flow | 2.69e-15 m³/s |
| Local stored-air balance | 6.63e-21 m³ |
| Global balance including ambient flow | 1.50e-19 m³ |

The tested closed tongue neck passes exactly zero flow. The mouth caps stay
clear of the complete solid-boundary sweep, with minimum clearance 5.55 mm
over these cases. During lip narrowing, two region pairs share one grid face;
whole-face events alone cannot certify each pair's opening. Persistent open
segments, untouched by any swept polygon boundary, independently certify
positive apertures for every internal pair throughout that interval. This
does not implement general contact scheduling for arbitrary shared-face
closures.

Accurate linear solves do not make a coarse time step accurate. The short lip
sections introduce fast undamped modes, and the coarser endpoint pressures
are not monotonically converged. Additional native runs refine the time step
beyond the full independent audit. Their final adjacent-refinement differences
are reported separately in the receipt; these finer runs are temporal
convergence evidence, not additional independent geometry audits.

| Interval | Final subdivision comparison over 1 ms | Maximum endpoint pressure difference |
|---|---|---:|
| Tongue closure | 512 → 1,024 | 0.002854 Pa |
| Tongue release | 512 → 1,024 | 0.007389 Pa |
| Lip narrowing | 1,024 → 2,048 | 0.006564 Pa |

These subdivisions resolve the declared undamped controls. They do not
establish the time resolution needed by a future complete lossy speech engine.

The frozen split/unsplit controls conserve energy within 8.89e-16 relative;
maximum pressure differences are 2.85e-14 Pa for the fixed tube and 1.82e-12 Pa
for the moving tube. The original two sealed tree controls reproduce their
historical states exactly. Complete-lip-closure rejection is also checked.

Release and AddressSanitizer/UndefinedBehaviorSanitizer builds pass the
coupling controls. The sanitizer run independently checks all three moving
cases at 16 subdivisions, the symmetry/closure controls, the original tree
regressions, and the geometry export/scaling controls. Its receipt is
`out/articulated-lip-asan/receipt.json`. Four focused native CTests pass:
core engine, planar airspace, moving face geometry and general-graph acoustics.
All 273 Release and 207 sanitizer recorded source, binary, input and artifact
hashes were verified after the completed runs.

## Reproduce and continue

```sh
cmake --build build --target ducklips_airspace_geometry_probe ducklips_planar_inertive_motion_probe
mkdir -p out/articulated-lip/poses
build/ducklips_airspace_geometry_probe out/nato-airspace-geometry.csv 5.020 5.085 512 out/articulated-lip/geometry.csv --frames out/articulated-lip/poses --core bodies/reference.body
python tools/planar_lip_reference.py
```

The geometry and airway CSV inputs are the saved native NATO replay; their
hashes and the input body, sources, binaries and generated fixtures are
recorded in the receipt. Historical receipts retain their original hashes.

The next speech integration work must handle disappearing lip storage,
connect the actual epilarynx and dynamic vocal folds at the inlet seam, and
restore the physical loss/radiation behavior at these new connections. The
coarse planar geometry also needs continuous state transport through all
region splits and merges. None of the pressure controls is listening evidence
or a demonstration of the requested rendering speed; performance optimization
remains parked.

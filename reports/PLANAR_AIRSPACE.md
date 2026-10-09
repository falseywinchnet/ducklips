# Air volume and openings from solid geometry

The [spatial-pressure speech experiment](SPATIAL_PRESSURE.md) exposed finite
air-volume jumps when the tongue tip crossed a section's measurement ray.
The user then [confirmed](LISTENING_SPATIAL_2026_10_09.md) that November was
improved in that recording, while Foxtrot remained jammed around /tr/.

This change adds native finite-volume geometry that accounts for air above
and below solid tongue parts. It is a geometric reference for the physical
engine, **not yet connected to speech rendering**. No newer audio or listening
score is attributed to it. The research goal remains realistic articulated
speech, including moving pressure/acoustic state and the eventual speed gate.

## What is measured

`PlanarAirspace` accepts a polygonal tract boundary, overlapping solid parts,
a grid and an extrusion depth, all in SI units. It measures:

```
cell volume = depth × area(cell ∩ (tract domain minus union of solids))
shared aperture = depth × length(shared face occupied by free air)
```

The union prevents body, blade and pad overlaps being subtracted twice.
The method retains all free intervals on a slice, including those beneath
the tongue. Its area calculation splits at polygon vertices, edge crossings,
and crossings of grid lines. Between these events the interval bounds are
linear, so midpoint integration is exact for the supplied polygons. This is
geometric integration, not a spatial raster average or a temporal filter.
Curved tissue is still polygonized and needs contour-refinement checks.

Preparation allocates bounded scratch space; builds allocate nothing. Internal
vertical openings require air to reach the face from both sides, so a solid
edge coincident with a grid face seals it. Outer grid faces are left zero:
wall/port boundary conditions must be specified by an eventual solver.

Cell volume and connecting aperture are separate quantities. A cell may
contain two disconnected pieces of air: adding their volumes does not make
them one pressure compartment. The API documents that restriction. The
independent audit checks for split cells and compares graph connectivity
with exact polygon connectivity before accepting its selected reference poses.
General component reconstruction is still required before arbitrary moving
anatomy can use one pressure unknown per geometrical cell. The subsequent
[66-pose transport audit](CAVITY_TRANSPORT.md) finds eight poses with split
cells, including one that would falsely join the sealed /t/ chambers.

## Controlled physical checks

The native test checks overlapping solids against independently known union
areas, crossing edges, a concave tract, reversed contour winding, refinement
of the grid, and barriers placed exactly on vertical or horizontal faces.
These barriers preserve both chambers' volumes and have exactly zero open
aperture across contact.

A rounded pad translates inside a closed rectangular mouth. All air above
and below it remains accounted for: total air volume agrees with domain area
minus pad area to 6.98e−20 m³. Reducing displacement tenfold reduces the maximum
individual cell-volume change tenfold; the old upper-ray jump did not vanish
with displacement. Circle polygons of 32, 64, 128 and 256 vertices converge
to the analytic circular area, with approximately fourfold error reduction
per doubling. These are checks of geometry and conservation, not human anatomy.

The full Release suite passed 14/14. The new native geometry test passed
ASan/UBSan, 1/1, with the existing `-O1 -g` sanitizer configuration. The
previous pressure/stream sanitizer results are not counted as fresh tests here.

## Actual Foxtrot poses

The native replay tool now also exports tongue-root posture. All existing
fields in all 23,239 rows remain identical to the geometry from the assessed
spatial-pressure recording. The geometry probe reconstructs 66 poses from
5.020–5.085 s, including the /t/→/r/ movement; maximum reconstructed attachment
and tip-coordinate error is 9.80e−12 cm, consistent with the saved CSV precision.

The solid reference uses the actual conservative body material map and pad
positions, the existing tract wall/floor profile, and an explicit 2.2 cm
extrusion. Its blade is a diagnostic solidification of the former line:
the upper attachment is the existing body material point; a lower attachment
starts 1.2 cm below it in the resting body and follows the same material map.
The blade tapers toward the rounded pad. This gives it a finite connection
to the body. Merely offsetting thickness outward from the upper attachment
left a point connection and a false air route through the attachment; that
construction was rejected before accepting the reference.

This extruded geometry replaces the empirical upper-gap power law and adds
an explicit blade volume. It therefore is **not an isolated numerical patch**
to the old anatomical model. The two geometries must not be presented as
acoustically equivalent or as independently measured anatomy.

| Control transition | Old upper-gap oral-volume change | Solid-airspace reference, 1,024 contour vertices |
|---|---:|---:|
| 5.024→5.025 s | +2.291309 mL | +0.021850 mL |
| 5.070→5.071 s | +2.731570 mL | +0.069394 mL |
| 5.072→5.073 s | +2.844607 mL | +0.056268 mL |

Across the 66 poses, maximum absolute volume difference from the 1,024-vertex
reference is 0.007880 mL at 128 vertices, 0.001886 mL at 256, and 0.000378 mL
at 512. The largest one-millisecond change in the finest reference is
0.083452 mL. Small geometric changes no longer produce the old multi-millilitre
sampling jumps in these tested poses. This does not yet establish an audible
improvement or validate the remaining volume changes physically.

## Independent volume, aperture and connectivity audit

Shapely/GEOS independently constructs the Boolean fluid domain for four
512-vertex snapshots at 5.025, 5.040, 5.071 and 5.073 s. No invalid polygons
are repaired automatically. It clips every one of the 4,096 cells and every
internal face independently of the native strip integration.

- Maximum cell-volume discrepancy: 3.54e−21 m³.
- Maximum total-volume discrepancy: 2.44e−19 m³.
- Maximum aperture-area discrepancy: 6.42e−18 m².
- No cell contains disconnected pieces of air in these four snapshots.
- Graph components using strictly positive native volumes/apertures match
  the exact polygon components and their volumes; no closure threshold is
  introduced to obtain that match.

At 5.040 s, the /t/ closure separates rear air (56.195761 mL) from front and
sublingual air (12.670920 mL). Before the closure and during the measured
/r/ approach, the air reconnects. The comparison preserves the sealed stop
while retaining air below the tongue. The full-domain diagram is saved as
`out/planar-airspace/foxtrot-closure.png`.

Commands, source/binary/input hashes, contour refinement, independent errors,
and the figure hash are in [planar_airspace_reference.json](planar_airspace_reference.json).

```sh
cmake --build build
ctest --test-dir build --output-on-failure
# Generate the full native geometry replay with the same alphabet text used
# in spatial_pressure_reference.json, retaining the new posture_root column.
python tools/planar_airspace_reference.py
```

The audit uses locally installed Shapely only as an independent offline
reference. The native geometry component adds no external C++ dependency.

## Remaining engine work

The current speech engine still uses the old section geometry. The
[mean-pressure solver now supports cyclic and disconnected graphs](CAVITY_GRAPH.md),
with fixed-pose pressure checks on these snapshots; the acoustic solver still
requires a tree.
The next integration must preserve actual fluid components, assign physical
ports and wall boundaries, derive inertance/resistance from the same geometry,
and transport stored mass/acoustic state conservatively as those components
move and close. The nasal path, existing side cavities and glottal/radiation
coupling must remain represented; this planar reference currently describes
the oral/pharyngeal domain only. Adding the lower volume to the old common pressure would
reintroduce the topological defect.

The reference's blade is not a constant-volume muscle model. Body-plus-blade
union volume is not constrained; the 2.2 cm extrusion is a stated modeling
choice rather than matched transverse anatomy. Contour approximations,
arbitrary split cells, lateral channels, physical tissue contact, moving
boundary work and boundary-condition integration remain unresolved. Those
requirements must survive integration; the small geometry errors here do
not prove realistic speech, a matched individual body, or 10× realtime.

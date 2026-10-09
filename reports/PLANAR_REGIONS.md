# Separate air regions inside a grid cell

The native planar geometry now resolves disconnected air pieces within each
cell and attaches openings to the correct pieces. All 66 saved Foxtrot poses
from 5.020–5.085 s agree with an independent polygon calculation, including
all eight poses that failed the previous one-pressure-per-cell audit.

This removes the spatial false connection at the /t/ seal in the reference
geometry. It does not yet supply moving-region identity, consistent geometry
over a time interval, or acoustic speech integration. There is no new audio
or listening score. [November improved; Foxtrot remained jammed around /tr/](LISTENING_SPATIAL_2026_10_09.md)
remains the latest user assessment of the spatial-pressure recording.

## Representation

`PlanarAirspace::build_regions` retains the existing ordinary cell-volume
and grid-face totals, and additionally emits:

- One `PlanarRegion` per connected air piece within a cell, with its volume
  and an interior representative point. The point is not a volume centroid.
- One `PlanarRegionLink` per region pair and grid face, containing the total
  aperture area. A grid face can carry links between several distinct region
  pairs. Its face index follows the existing vertical/horizontal face arrays.

The geometry is divided into vertical strips at polygon vertices,
intersections and grid crossings. Within each strip, free-space bounds are
linear. The region calculation retains the polygon-edge identity of those
bounds so it can evaluate both ends of a strip. Positive-volume pieces in the
same cell join only through a shared interval of positive length. Contact at
one point does not connect air regions. A grid boundary produces a port
instead of joining the two regions into a common pressure compartment.

Horizontal openings use the identity of the common free interval. This
avoids trying to infer the connection from rounded endpoint equalities.
After within-cell unions, ports belonging to the same pair of regions on the
same face are summed. No aperture floor or small-fragment deletion is used.

`prepare_regions` allocates bounded scratch capacity. Building regions does
not allocate; insufficient capacity returns failure. Preparing a different
grid invalidates the old region workspace, and an ordinary geometry build
clears stale region outputs. Region indices are local to the current pose:
reusing an index as a moving pressure-state identity is not supported.

The existing solid-union rules, transverse-depth assumption, and unvalidated
blade model are unchanged. The output describes the stated planar reference,
not independently measured three-dimensional human anatomy.

## Independent verification of every region and opening

The probe accepts the exact saved contours and grid. Shapely/GEOS independently
constructs the Boolean fluid domain and clips it to every cell. Each native
region must match exactly one clipped fluid component through its interior
point. Counts must agree in every cell; matching total cell volume alone is
insufficient. Every expected opening between neighbouring components is then
constructed geometrically, including its identity and area.

The audit also checks that each connected native graph component belongs to
exactly one whole-domain fluid component. Matching component counts and total
volumes alone would not detect all erroneous exchanges of pieces.

Across all 66 poses:

| Check | Maximum discrepancy |
|---|---:|
| Individual region volume | 4.49e−21 m³ |
| Individual opening area | 2.76e−18 m² |
| Whole connected air-region volume | 3.94e−19 m³ |

Every component count, component assignment and opening identity agrees.
There are 1,269–1,298 native regions per pose. All eight previously rejected
split-cell poses are resolved: 5.032, 5.037 and 5.079–5.084 s. These results
cover the supplied poses and controlled tests, not arbitrary future anatomy
or contacts below floating-point resolution.

## The /t/ seal in pressure

At 5.032 s, one grid cell contains fragments from opposite sides of the stop.
The old combined-cell graph had one connected air region; the region graph
has the correct two. A pressure control initializes the rear chamber at
600 Pa and the front chamber at zero, with a sealed outer boundary.

Both representations use the same explicit linear-mobility fixture,
1e−7 m²/(Pa s), and a 1 ms pressure step. This mobility is a numerical control,
not calibrated vocal-tract resistance.

| Representation | Result |
|---|---|
| Separate regions | Zero pressure change in both chambers; every port flow is exactly zero. |
| Combined grid cells | Initial mixing transfers 1.99e−14 m³ of reference-volume storage into the front chamber; after the step the gain is 6.32e−12 m³ and mean front pressure is 0.07325 Pa. |

The combined model conserves its total stored amount but transfers it across
a sealed boundary. Thus global conservation alone did not establish correct
contact behavior. The region model preserves both conservation and isolation
in this fixed-pose check. This comparison does not quantify an audible change
in the production speech engine, which still uses section geometry.

## Native tests and reproduction

Controlled native tests cover an internal wall dividing one cell, two pieces
of tissue meeting at a single point, a 1e−12 m thick wall, a 1e−10 m opening,
separate ports sharing either orientation of grid face, explicit capacity
failure, repeated preparation, stale-output invalidation, and allocation-free
component building. Earlier volume and aperture tests remain in place.

The full Release suite passed 14/14. The final workspace-lifecycle changes
then passed focused Release and ASan/UBSan geometry tests. The final geometry
test passed 1/1 under sanitizers. All 66 native geometry probes were also run
under ASan/UBSan and independently checked against GEOS; every case passed.
That sanitizer sweep concerns geometry; its unchanged pressure-control
helper uses the Release pressure probe.

[planar_regions_reference.json](planar_regions_reference.json) records every
pose, command, input/output hash, source/binary hash, dependency version and
the pressure control. Full sanitizer evidence is retained locally in
`out/planar-regions-asan/receipt.json`.

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/planar_regions_reference.py
python tools/planar_regions_reference.py \
  --binary .build/asan/ducklips_planar_regions_probe \
  --output out/planar-regions-asan \
  --report out/planar-regions-asan/receipt.json
```

The audit consumes the 66 snapshots exported by the command recorded in
[cavity_transport_reference.json](cavity_transport_reference.json).
The native implementation adds no external C++ geometry dependency.

## Next integration requirement

The previous [transport experiments](CAVITY_TRANSPORT.md) still have unresolved
temporal refinement. A correct spatial decomposition does not establish how
old air content, pressure or momentum maps into regions that split, merge,
appear and disappear. Assigning state by the new output index would again
transfer air incorrectly. Volumes, openings and state correspondence must
come from consistent geometry over the step.

After that, physical mean-flow laws, acoustic state and moving-boundary work
must preserve the nasal path, side cavities and glottal/radiation coupling.
Tissue-volume conservation, transverse anatomy and new listening evidence
also remain required. The realism and eventual speed gates remain unchanged;
performance optimization is still parked.

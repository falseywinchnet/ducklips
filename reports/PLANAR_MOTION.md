# Measure the moving geometry, not blended endpoint metrics

The native planar geometry can now sample a declared motion between two
corresponding contour sets and recompute regions, volumes and openings from
that geometry. It no longer needs to substitute interpolated endpoint
measurements for intermediate geometry in this research path.

The complete 65 ms Foxtrot interval was checked at 125 µs spacing: 521 poses.
Every sampled region, opening, and whole-fluid component agrees with the
independent polygon calculation. This includes intermediate split cells that
were absent at both saved millisecond endpoints.

This is a geometry advance. Conservative pressure-state correspondence across
moving regions and acoustic integration remain unfinished; there is no new
speech recording or listening score. The latest assessment remains
[November improved, Foxtrot jammed around /tr/](LISTENING_SPATIAL_2026_10_09.md).

## Native operation and its assumptions

`PlanarAirspace::build_between` receives two matched contour sets and a fraction
in [0,1]. It moves corresponding vertices along straight paths, then calls the
same component-resolved geometric extraction used for the endpoint poses.
It does not interpolate cell volume, face aperture, or component indices.
Preparation allocates the sampled-contour workspace; sampling allocates nothing.
Nonfinite/out-of-range fractions, mismatched contour sizes and insufficient
capacity fail explicitly. The sampled contours are available for diagnostics.

The caller supplies vertex correspondence and valid intermediate polygons.
Endpoint storage must remain independent of the object's sampled workspace.
For this experiment the saved domain, body, blade and pad vertices retain
their original ordering. The reference uses the existing 2.2 cm extrusion.
The piecewise linear vertex path is an explicit assumption: it does not
reconstruct the exact articulator dynamics, enforce constant tissue volume,
or establish a muscle model. Those anatomical requirements remain open.

## A coupled analytical check

A sloping solid boundary sweeps a triangular air wedge into a cell. With
height h, its volume is h²/2 and its opening into the adjacent region is h
(unit extrusion depth). At the midpoint of h=0→0.4, the actual cell volume is
0.02; interpolation of endpoint volumes gives 0.04.

The test uses the native sampled volume and aperture in `CavityNetwork`, with
the adjacent region represented by a zero-pressure reservoir. For this
linear height motion, midpoint conductance is the exact interval average.
The independently derived finite pressure is recovered through both creation
and disappearance at 1, 2, 4 and 8 timesteps. Maximum pressure error is
2.98e−14 Pa and reference-volume balance error is 6.81e−16 m³. These are
controlled numerical geometries, not measurements of a human mouth.

This demonstrates why a conserving pressure equation is insufficient when
its volume and aperture histories describe different motions.

## What endpoint blending missed in Foxtrot

The 521 samples use the exact saved poses as interval endpoints and eight
subintervals per millisecond. The native sampled vertex positions agree with
an independent linear interpolation to 1.39e−17 m. Shapely/GEOS then audits
the actual native contours, not a separately rounded copy of the path.
All polygons are valid and none are automatically repaired.

Maximum discrepancies against the independent geometry are:

| Quantity | Maximum discrepancy |
|---|---:|
| Individual region volume | 4.50e−21 m³ |
| Region-pair opening area | 4.28e−18 m² |
| Whole fluid-component volume | 4.07e−19 m³ |

The comparison with blended endpoint metrics exposes substantial local errors:

| Saved-pose time convention | Moving-contour measurement | Endpoint-metric interpolation |
|---|---:|---:|
| 5.0675 s, cell 3688 | 0.001375 mL | 0.013053 mL |
| 5.07025 s, internal vertical face 3743 | Exactly closed | 0.270703 cm² open |

That false opening is 75% of the entire grid-face area. It is an internal
computational face, not a measured anatomical constriction. The discrepancy
establishes that metric interpolation can invent a passage in the declared
geometry. Its audible contribution has not been measured.

There are 78 sampled poses with split cells. Sixteen intermediate samples,
spread over 12 original millisecond intervals, contain split cells even though
both interval endpoints have none. Examples include 5.028625, 5.030125,
5.030250, 5.049500 and 5.049625 s. Thus endpoint-only acceptance of a
one-pressure-per-cell representation is insufficient, including in the
previous partial transport experiments.

Whole-domain closure/reopening is bracketed by the sampled geometry:

- One chamber at 5.031250 s; two by 5.031375 s.
- Two chambers at 5.062375 s; one by 5.062500 s.

These brackets describe the declared interpolated geometry and the saved
pose time labels. They are not direct measurements of the speaker or a new
rendered speech waveform. A finite sample grid also does not prove the absence
of shorter events between samples.

## Validation and evidence

The full Release suite passed 14/14. The geometry test, including the new
coupled motion/pressure checks, passed ASan/UBSan, 1/1. Eleven native motion
cases were additionally run under sanitizers and independently checked against
GEOS: endpoints, hidden splits, both contact brackets, and the two largest
metric-blending discrepancies. All passed. This turn did not rerun streaming
under sanitizers; the fresh sanitizer scope is geometry and its coupled fixture.

[planar_motion_reference.json](planar_motion_reference.json) records all 521
samples, exact commands, input/output hashes, native source/binary hashes,
reference-library versions and discrepancy extrema. Each sampled geometry also
has saved cell/face metrics for subsequent integration work. Sanitizer receipts
are in `out/planar-motion-asan/receipt.json`.

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/planar_motion_reference.py
python tools/planar_motion_sanitizer.py
```

The shared independent geometry checker was extracted into `check_geometry`
in `tools/planar_regions_reference.py`, so static and moving audits use the
same complete region/port verification. No external C++ dependency was added.

## Remaining integration requirement

The next pressure integration must account for openings over each time interval
and carry stored air consistently through region continuation, splitting,
merging, appearance and disappearance. A sampled output index does not identify
the same air region at the next instant. Neither arbitrary index reuse nor
independent interpolation of geometric metrics is justified by these checks.

The previous transport model's pressure-refinement failure remains unresolved.
Physical flow laws, acoustic state and moving-boundary work must still retain
the nasal path, side cavities and glottal/radiation coupling. Tissue-volume
conservation, individual transverse anatomy and new listening evidence remain
requirements. Performance work remains parked and no realism or speed gate
has been declared satisfied.

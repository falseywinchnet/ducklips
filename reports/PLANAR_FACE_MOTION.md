# Preserve opening and closure times within a motion step

`PlanarFaceMotion` now measures a fixed grid face over the declared linear
motion of contour vertices. It returns chronological intervals and the
integrated opening in each interval. The output retains intervening closures,
including a face that is open at both saved endpoints.

This supplies temporal geometry for the pending moving-pressure integration.
It does not yet carry pressure between moving regions or render speech.
Foxtrot's audible /tr/ problem remains unresolved; no new recording or
listening score is claimed.

## Why chronology matters

The previous motion audit found that Foxtrot's internal vertical face 3743
was closed at 5.07025 s, while endpoint interpolation left it 75% open.
The new calculation follows this face throughout 5.070–5.071 s:

| Quantity | Declared moving contours | Endpoint-metric approximation |
|---|---:|---:|
| Closure time, saved-pose convention | 5.0701887716 s | Only reaches zero at 5.071 s |
| Fraction of the interval fully closed | 81.1228% | 0%, except the final instant |
| Average face opening | 0.0450169 cm² | 0.180469 cm² |

This is an internal computational face in the existing 2.2 cm extrusion,
not a measured anatomical constriction. These times describe the declared
linear vertex path and saved geometry labels, not measured human contact
or an observed acoustic event. They should not be read as anatomical timing
precision.

A coupled pressure control demonstrates why retaining only an accurate mean
area is still insufficient. Two fixed chambers start at zero pressure. Their
connecting shutter closes halfway through a step; a source supplies air to
one chamber afterward. Processing the native opening windows in order gives
exactly zero transfer to the other chamber and preserves the supplied air.
Using one whole-step mean opening instead produces a false 0.28185 Pa pressure
rise in the other chamber. The test uses an explicit conductance proportional
to area, not a calibrated human-flow law.

## Native method

For one face, project a moving polygon edge into normal and tangent
coordinates. Its intersection with the supporting line has tangent coordinate
N(t)/D(t), where N has degree at most two and D at most one. The ordering and
clipping of these intersections can change at:

- A vertex crossing the face's supporting line: a linear root.
- An edge intersection reaching a face endpoint: a quadratic root.
- Two edge intersections exchanging order: a root of
  N₁D₂ − N₂D₁, of degree at most three.
- A denominator becoming zero, which is also retained as a cut.

The implementation partitions normalized time at these candidates. Edge
bounds discard candidates that cannot reach the face; algebraic roots may
still add harmless extra cuts outside the actual segments. Nearby roots are
not merged using a distance tolerance. Roots are isolated between derivative
roots, with near-zero stationary points retained for tangency cases.

Within each interval, intersections from both sides of the face determine
the domain opening minus the union of solids. This preserves a seal when a
solid boundary lies exactly on the face. Adaptive 4/8 point Gauss quadrature
integrates the resulting aperture between the event cuts. The returned
quadrature discrepancy is an estimate, not a rigorous error bound. Polynomial
root finding and geometry remain floating-point calculations, not certified
exact predicates.

Preparation reserves bounded workspace. The build operation allocates no
storage. Invalid inputs, exhausted event/edge capacity, and quadrature that
cannot meet its target fail explicitly; all outputs must then be discarded.
The caller supplies matching vertex identities and simple intermediate
polygons. The method neither infers body motion nor enforces tissue-volume
conservation.

Each output area sums the openings across one geometric face. It does not
identify the region pair on either side. Disconnected openings must not be
silently coupled as one pressure passage. Nor does average area determine
average resistance for a nonlinear physical flow law.

## Verification

Native analytical controls cover:

- A face open at both endpoints but closed for one third of the interval.
- A closure lasting only 10⁻⁷ of the normalized interval.
- A rational aperture with an independently integrated logarithmic solution.
- Two overlapping solids exchanging ownership of the open boundary.
- A cubic ordering equation with roots at 0.2, 0.5 and 0.8, retaining both
  separate open intervals. Integral discrepancy: 3.97e−15 m².
- Time reversal, rotation, winding independence and coincident boundaries.
- The post-closure pressure-source control, conservation, bounded capacities
  and no allocation during evaluation.

The independent audit uses Shapely/GEOS to measure the actual face opening
and SciPy to integrate it in time. Its quadrature receives uniform time
partitions, never native event roots. It refines the independent partition
when QUADPACK reports roundoff trouble; warnings and intermediate estimates
are retained in the receipt. The initial relative tolerance was also found
too loose for the larger synthetic geometries and tightened to 1e−12 before
the final audit. No native tolerance or geometry was changed to fit the oracle.

All 49 cases passed: 25 Foxtrot face intervals and 24 pairs of independently
moving convex triangles. Maximum native/oracle integral discrepancy was
1.20e−12 m²; maximum difference between the final two independent refinements
was 1.80e−13 m². Three Foxtrot cases needed further independent subdivision
after a warning, up to 2048 uniform subintervals. Every final independent
integration completed without a warning. Three interior probes of every
reported closed run found zero opening within 1e−16 m².

The full Release suite passed 15/15. The planar geometry and new face-motion
tests passed ASan/UBSan, 2/2. A separate sanitizer audit passed seven cases
(five Foxtrot intervals and two synthetic motions), each with the independent
GEOS/SciPy checks; maximum discrepancy was 6.28e−13 m². Streaming was not
rerun under sanitizers in this change. All recorded source, binary, input and
output hashes were checked against the final files.

These checks cover the selected faces and analytical controls. They do not
constitute an audit of every face in every interval, a proof of exact root
predicates, or a moving pressure-state integration.

[planar_face_motion_reference.json](planar_face_motion_reference.json) records
all cases, commands, input/output hashes, source/binary hashes, library versions,
independent refinement results and closed intervals. Full native windows are
saved in `out/planar-face-motion/`. The sanitizer audit has its own receipt in
`out/planar-face-motion-asan/receipt.json`.

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/planar_face_motion_reference.py
# A single saved interval, using the existing pair-input format:
./build/ducklips_planar_regions_probe --face-v 3743 \
  < out/planar-motion/interval-050.txt
```

## Remaining speech integration

This resolves a temporal-geometry requirement, not the complete moving-airway
problem. Region continuation, splitting/merging, disappearance, stored-air
transfer, physical flux laws, acoustic state and moving-wall work remain to
be connected. Opening histories must be consumed in chronological order with
those states; treating their average as a universal passage parameter would
reintroduce the demonstrated error.

The production spatial-pressure engine still uses its existing section
geometry and tree acoustics. The nasal path, side cavities, glottal coupling
and radiation must be preserved by the eventual replacement. Transverse
anatomy, tissue-volume conservation and new human listening evidence remain
unverified. Performance work remains parked.
